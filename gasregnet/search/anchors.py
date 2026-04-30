"""Profile-driven anchor detection over RefSeq catalogs."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import duckdb
import polars as pl

from gasregnet.config import GasRegNetConfig, load_config
from gasregnet.datasets.refseq import read_refseq_catalog_manifest
from gasregnet.schemas import AnchorHitsSchema, validate
from gasregnet.search.hmmer import hmmsearch

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


def _empty_anchor_hits() -> pl.DataFrame:
    return validate(pl.DataFrame(schema=ANCHOR_HITS_SCHEMA), AnchorHitsSchema)


def _profile_path(profile_dir: Path, anchor_family: str) -> Path:
    return profile_dir / f"{anchor_family}.hmm"


def _feature_metadata(db: Path, protein_accessions: list[str]) -> pl.DataFrame:
    if not protein_accessions:
        return pl.DataFrame(
            schema={
                "protein_accession": pl.Utf8,
                "locus_tag": pl.Utf8,
                "gene": pl.Utf8,
                "product": pl.Utf8,
            },
        )
    with duckdb.connect(str(db), read_only=True) as connection:
        return connection.execute(
            """
            select
                proteins.protein_accession,
                coalesce(features.locus_tag, '') as locus_tag,
                coalesce(features.gene, '') as gene,
                coalesce(features.product, proteins.description) as product
            from proteins
            left join features using (protein_accession)
            where proteins.protein_accession in (
                select unnest(?::varchar[])
            )
            """,
            [protein_accessions],
        ).pl()


def _hits_for_family(
    *,
    dataset_name: str,
    db: Path,
    protein_faa: Path,
    analyte: str,
    anchor_family: str,
    profile_hmm: Path,
    e_value_threshold: float,
    bitscore_threshold: float | None,
) -> pl.DataFrame:
    hits = hmmsearch(profile_hmm, protein_faa, e_value=e_value_threshold)
    if hits.is_empty():
        return _empty_anchor_hits()
    hits = hits.filter(pl.col("evalue") <= e_value_threshold)
    if bitscore_threshold is not None:
        hits = hits.filter(pl.col("score") >= bitscore_threshold)
    if hits.is_empty():
        return _empty_anchor_hits()

    metadata = _feature_metadata(db, hits["target_id"].unique().to_list())
    frame = (
        hits.join(
            metadata,
            left_on="target_id",
            right_on="protein_accession",
            how="left",
        )
        .with_columns(
            [
                pl.lit(dataset_name).alias("dataset_name"),
                pl.lit(analyte).alias("analyte"),
                pl.lit(anchor_family).alias("anchor_family"),
                pl.col("target_id").alias("protein_accession"),
                pl.col("locus_tag").fill_null("").alias("locus_tag"),
                pl.col("gene").fill_null("").alias("gene"),
                pl.col("product").fill_null("").alias("product"),
                pl.col("score").alias("bitscore"),
                pl.col("evalue").alias("e_value"),
                pl.lit(None, dtype=pl.Float64).alias("identity"),
                pl.lit(None, dtype=pl.Float64).alias("coverage"),
                pl.lit("hmmer").alias("evidence_type"),
            ],
        )
        .select(list(ANCHOR_HITS_SCHEMA))
        .unique(
            subset=[
                "dataset_name",
                "analyte",
                "anchor_family",
                "protein_accession",
                "locus_tag",
            ],
        )
    )
    return validate(frame, AnchorHitsSchema)


def detect_anchors_profile(
    manifest: Path,
    *,
    config: Path | GasRegNetConfig,
    root: Path,
    profile_dir: Path,
    bitscore_threshold: float | None = None,
    e_value_threshold: float = 1e-20,
) -> pl.DataFrame:
    """Run HMM profile anchor detection against RefSeq catalog protein FASTAs."""

    catalogs = read_refseq_catalog_manifest(manifest, root=root)
    loaded_config = load_config(config) if isinstance(config, Path) else config
    frames: list[pl.DataFrame] = []
    for catalog in catalogs.iter_rows(named=True):
        dataset_name = str(catalog["dataset_name"])
        db = Path(str(catalog["out_db"]))
        protein_faa = Path(str(catalog["protein_faa"]))
        for analyte in loaded_config.analytes:
            for family in analyte.anchor_families:
                profile_hmm = _profile_path(profile_dir, family.name)
                family_hits = _hits_for_family(
                    dataset_name=dataset_name,
                    db=db,
                    protein_faa=protein_faa,
                    analyte=analyte.analyte,
                    anchor_family=family.name,
                    profile_hmm=profile_hmm,
                    e_value_threshold=e_value_threshold,
                    bitscore_threshold=bitscore_threshold,
                )
                if not family_hits.is_empty():
                    frames.append(family_hits)
    if not frames:
        return _empty_anchor_hits()
    return validate(pl.concat(frames, how="vertical"), AnchorHitsSchema)
