"""DIAMOND-backed seed rescue for sharded corpus-store searches."""

from __future__ import annotations

import logging
import tempfile
from collections.abc import Iterable
from pathlib import Path
from typing import Any, Literal, Protocol, cast

import duckdb
import polars as pl

from gasregnet.config import AnalyteConfig
from gasregnet.datasets.corpus_store import register_views
from gasregnet.errors import MissingInputError
from gasregnet.schemas import AnchorHitsSchema, validate
from gasregnet.search import diamond

LOGGER = logging.getLogger(__name__)

ANCHOR_HITS_SCHEMA: dict[str, Any] = {
    "dataset_name": pl.Utf8,
    "analyte": pl.Utf8,
    "anchor_family": pl.Utf8,
    "protein_accession": pl.Utf8,
    "locus_tag": pl.Utf8,
    "gene": pl.Utf8,
    "product": pl.Utf8,
    "bitscore": pl.Float64,
    "e_value": pl.Float64,
    "identity": pl.Float64,
    "coverage": pl.Float64,
    "evidence_type": pl.Utf8,
}


class ShardLike(Protocol):
    """Structural shard type accepted by seed rescue."""


class CorpusStoreLike(Protocol):
    """Structural corpus store handle accepted by seed rescue."""


def _empty_anchor_hits() -> pl.DataFrame:
    return validate(pl.DataFrame(schema=ANCHOR_HITS_SCHEMA), AnchorHitsSchema)


def _attr(obj: object, *names: str) -> object | None:
    for name in names:
        if hasattr(obj, name):
            return cast(object, getattr(obj, name))
    if isinstance(obj, dict):
        for name in names:
            if name in obj:
                return cast(object, obj[name])
    return None


def _string_list(value: object | None) -> list[str]:
    if value is None:
        return []
    if isinstance(value, str):
        return [value]
    if isinstance(value, Iterable):
        return [str(item) for item in value]
    return [str(value)]


def _frame_from_store_method(store: object, shard: object) -> pl.DataFrame | None:
    for method_name in (
        "fetch_proteins_for_shard",
        "fetch_shard_proteins",
        "proteins_for_shard",
    ):
        method = getattr(store, method_name, None)
        if callable(method):
            return cast(pl.DataFrame, method(shard))

    fetch_dataset = getattr(store, "fetch_proteins_for_dataset", None)
    if callable(fetch_dataset):
        frames = []
        phylum = _attr(shard, "phylum")
        for dataset_name in _string_list(
            _attr(shard, "dataset_names", "datasets", "dataset_name"),
        ):
            try:
                frames.append(fetch_dataset(dataset_name, phylum))
            except TypeError:
                frames.append(fetch_dataset(dataset_name))
        if frames:
            return cast(pl.DataFrame, pl.concat(frames, how="diagonal_relaxed"))
    return None


def _frame_from_partitioned_store(store_root: Path, shard: object) -> pl.DataFrame:
    dataset_names = _string_list(
        _attr(shard, "dataset_names", "datasets", "dataset_name"),
    )
    phylum = _attr(shard, "phylum")
    predicates: list[str] = []
    params: list[object] = []
    if dataset_names:
        predicates.append("proteins.dataset_name in (select unnest(?::varchar[]))")
        params.append(dataset_names)
    if phylum is not None:
        predicates.append("proteins.phylum = ?")
        params.append(str(phylum))
    where_clause = f"where {' and '.join(predicates)}" if predicates else ""
    with duckdb.connect() as connection:
        register_views(connection, store_root)
        return connection.execute(
            f"""
            select
                proteins.dataset_name,
                proteins.protein_accession,
                coalesce(features.locus_tag, '') as locus_tag,
                coalesce(features.gene, '') as gene,
                coalesce(features.product, proteins.description, '') as product,
                proteins.sequence
            from proteins
            left join features using (dataset_name, protein_accession)
            {where_clause}
            """,
            params,
        ).pl()


def _proteins_for_shard(
    store: CorpusStoreLike | Path,
    shard: ShardLike,
) -> pl.DataFrame:
    proteins: pl.DataFrame
    if isinstance(store, Path):
        proteins = _frame_from_partitioned_store(store, shard)
    else:
        fetched = _frame_from_store_method(store, shard)
        if fetched is None:
            store_root = _attr(store, "store_root", "root", "path")
            if store_root is None:
                raise TypeError(
                    "store must expose fetch_proteins_for_shard, "
                    "fetch_proteins_for_dataset, or a partitioned store path",
                )
            proteins = _frame_from_partitioned_store(Path(str(store_root)), shard)
        else:
            proteins = fetched

    if proteins.is_empty():
        return proteins
    required = {"protein_accession", "sequence"}
    missing = required.difference(proteins.columns)
    if missing:
        raise ValueError(f"protein shard frame missing column(s): {', '.join(missing)}")
    if "dataset_name" not in proteins.columns:
        dataset_name = _attr(shard, "dataset_name")
        proteins = proteins.with_columns(
            pl.lit(str(dataset_name or "")).alias("dataset_name"),
        )
    for column in ("locus_tag", "gene", "product"):
        if column not in proteins.columns:
            proteins = proteins.with_columns(pl.lit("").alias(column))
    return proteins.select(
        [
            "dataset_name",
            "protein_accession",
            "locus_tag",
            "gene",
            "product",
            "sequence",
        ],
    )


def _write_shard_fasta(proteins: pl.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for row in proteins.iter_rows(named=True):
            sequence = str(row["sequence"]).replace("*", "")
            if not sequence:
                continue
            accession = str(row["protein_accession"])
            handle.write(f">{accession}\n")
            for offset in range(0, len(sequence), 80):
                handle.write(f"{sequence[offset : offset + 80]}\n")


def _seed_db_path(
    seeds_diamond_dir: Path,
    *,
    analyte: str,
    anchor_family: str,
) -> Path:
    return seeds_diamond_dir / f"{analyte}__{anchor_family}.dmnd"


def _threshold_percent(value: float) -> float:
    return value * 100.0 if value <= 1.0 else value


def _threshold_fraction(value: float) -> float:
    return value / 100.0 if value > 1.0 else value


def _anchor_hits_from_diamond(
    *,
    hits: pl.DataFrame,
    proteins: pl.DataFrame,
    analyte: str,
    anchor_family: str,
    identity_threshold: float,
    coverage_threshold: float,
) -> pl.DataFrame:
    if hits.is_empty():
        return _empty_anchor_hits()

    identity_percent = _threshold_percent(identity_threshold)
    coverage_percent = _threshold_percent(coverage_threshold)
    filtered = (
        hits.filter(
            (pl.col("percent_identity") >= identity_percent)
            & (pl.col("qcovhsp") >= coverage_percent),
        )
        .sort(
            ["query_id", "bitscore", "evalue", "percent_identity", "qcovhsp"],
            descending=[False, True, False, True, True],
        )
        .unique(subset=["query_id"], keep="first", maintain_order=True)
    )
    if filtered.is_empty():
        return _empty_anchor_hits()

    metadata = proteins.unique(subset=["protein_accession"], keep="first")
    frame = (
        filtered.join(
            metadata,
            left_on="query_id",
            right_on="protein_accession",
            how="inner",
        )
        .with_columns(
            [
                pl.lit(analyte).alias("analyte"),
                pl.lit(anchor_family).alias("anchor_family"),
                pl.col("query_id").alias("protein_accession"),
                pl.col("locus_tag").fill_null("").alias("locus_tag"),
                pl.col("gene").fill_null("").alias("gene"),
                pl.col("product").fill_null("").alias("product"),
                pl.col("evalue").alias("e_value"),
                (pl.col("percent_identity") / 100.0).alias("identity"),
                (pl.col("qcovhsp") / 100.0).alias("coverage"),
                pl.lit("seed_back_confirmed").alias("evidence_type"),
            ],
        )
        .select(list(ANCHOR_HITS_SCHEMA))
    )
    return validate(frame, AnchorHitsSchema)


def seed_rescue_for_shard(
    *,
    shard: ShardLike,
    store: CorpusStoreLike | Path,
    seeds_diamond_dir: Path,
    analytes: list[AnalyteConfig],
    identity_threshold: float = 0.30,
    coverage_threshold: float = 0.50,
    e_value: float = 1e-10,
    threads: int = 8,
    work_dir: Path | None = None,
    missing_binary: Literal["skip", "raise"] = "skip",
) -> pl.DataFrame:
    """Run DIAMOND seed rescue for one shard and return AnchorHitsSchema rows.

    ``store`` is intentionally structural. It may be a partitioned corpus-store
    path or an object exposing ``fetch_proteins_for_shard`` /
    ``fetch_proteins_for_dataset``. This lets the Phase F shard/store types wire
    in without making this module own their concrete definitions.
    """

    proteins = _proteins_for_shard(store, shard)
    if proteins.is_empty():
        return _empty_anchor_hits()

    temp_context = (
        tempfile.TemporaryDirectory(prefix="gasregnet_seed_rescue_")
        if work_dir is None
        else None
    )
    base_dir = Path(temp_context.name) if temp_context is not None else work_dir
    assert base_dir is not None
    try:
        query_value = _attr(shard, "protein_faa", "fasta_path", "query_faa")
        query_faa = Path(str(query_value)) if query_value is not None else None
        if query_faa is None or not query_faa.exists():
            query_faa = base_dir / "shard_proteins.faa"
            _write_shard_fasta(proteins, query_faa)

        frames: list[pl.DataFrame] = []
        for analyte_config in analytes:
            analyte = analyte_config.analyte
            for family in analyte_config.anchor_families:
                db = _seed_db_path(
                    seeds_diamond_dir,
                    analyte=analyte,
                    anchor_family=family.name,
                )
                if not db.exists():
                    LOGGER.info("skipping seed rescue; missing DIAMOND db: %s", db)
                    continue
                out_tsv = base_dir / f"{analyte}__{family.name}.tsv"
                try:
                    diamond.run_diamond(
                        query_faa,
                        db,
                        out_tsv,
                        evalue=e_value,
                        coverage=_threshold_fraction(coverage_threshold),
                        identity=_threshold_fraction(identity_threshold),
                        threads=threads,
                        sensitivity="ultra-sensitive",
                    )
                except MissingInputError as error:
                    if "diamond binary" in str(error) and missing_binary == "skip":
                        LOGGER.warning("skipping seed rescue; %s", error)
                        return _empty_anchor_hits()
                    raise

                if out_tsv.stat().st_size == 0:
                    continue
                hits = diamond.parse_diamond_output(out_tsv)
                rescued = _anchor_hits_from_diamond(
                    hits=hits,
                    proteins=proteins,
                    analyte=analyte,
                    anchor_family=family.name,
                    identity_threshold=identity_threshold,
                    coverage_threshold=coverage_threshold,
                )
                if not rescued.is_empty():
                    frames.append(rescued)

        if not frames:
            return _empty_anchor_hits()
        return validate(
            pl.concat(frames, how="vertical").unique(
                subset=[
                    "dataset_name",
                    "analyte",
                    "anchor_family",
                    "protein_accession",
                    "locus_tag",
                ],
            ),
            AnchorHitsSchema,
        )
    finally:
        if temp_context is not None:
            temp_context.cleanup()
