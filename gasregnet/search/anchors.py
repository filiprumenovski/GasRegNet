"""Profile-driven anchor detection over RefSeq catalogs."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import duckdb
import polars as pl

from gasregnet.config import AnalyteConfig, GasRegNetConfig, load_config
from gasregnet.datasets.refseq import read_refseq_catalog_manifest
from gasregnet.io.fasta import read_fasta
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
                "sequence": pl.Utf8,
            },
        )
    with duckdb.connect(str(db), read_only=True) as connection:
        return connection.execute(
            """
            select
                proteins.protein_accession,
                coalesce(features.locus_tag, '') as locus_tag,
                coalesce(features.gene, '') as gene,
                coalesce(features.product, proteins.description) as product,
                proteins.sequence
            from proteins
            left join features using (protein_accession)
            where proteins.protein_accession in (
                select unnest(?::varchar[])
            )
            """,
            [protein_accessions],
        ).pl()


def _all_protein_metadata(db: Path) -> pl.DataFrame:
    with duckdb.connect(str(db), read_only=True) as connection:
        return connection.execute(
            """
            select
                proteins.protein_accession,
                coalesce(features.locus_tag, '') as locus_tag,
                coalesce(features.gene, '') as gene,
                coalesce(features.product, proteins.description) as product,
                proteins.sequence
            from proteins
            left join features using (protein_accession)
            """,
        ).pl()


def _gene_symbol(description: str) -> str:
    for token in description.split():
        if token.startswith("GN="):
            return token.removeprefix("GN=")
    return ""


def _seed_for_family(analyte: AnalyteConfig, anchor_family: str) -> str:
    records = list(read_fasta(analyte.anchor_seeds))
    family_lower = anchor_family.lower()
    for _accession, description, sequence in records:
        if _gene_symbol(description).lower() == family_lower:
            return sequence
    for _accession, description, sequence in records:
        if family_lower in description.lower():
            return sequence
    return records[0][2]


def _seed_sequences(config: GasRegNetConfig) -> dict[str, dict[str, str]]:
    seeds: dict[str, dict[str, str]] = {}
    for analyte in config.analytes:
        seeds[analyte.analyte] = {
            family.name: _seed_for_family(analyte, family.name)
            for family in analyte.anchor_families
        }
    return seeds


def _sequence_match(query: str, seed: str) -> tuple[float, float]:
    query = query.replace("*", "")
    seed = seed.replace("*", "")
    if not query or not seed:
        return 0.0, 0.0
    aligned = min(len(query), len(seed))
    matches = sum(1 for left, right in zip(query, seed, strict=False) if left == right)
    identity = matches / aligned if aligned else 0.0
    coverage = aligned / max(len(query), len(seed))
    return identity, coverage


def _back_confirm(
    *,
    sequence: object,
    analyte: str,
    anchor_family: str,
    seeds_by_analyte: dict[str, dict[str, str]],
    identity_threshold: float,
    coverage_threshold: float,
) -> dict[str, object]:
    best_family = ""
    best_identity = 0.0
    best_coverage = 0.0
    for seed_family, seed_sequence in seeds_by_analyte.get(analyte, {}).items():
        identity, coverage = _sequence_match(str(sequence), seed_sequence)
        if (identity, coverage) > (best_identity, best_coverage):
            best_family = seed_family
            best_identity = identity
            best_coverage = coverage
    evidence_type = "profile_match"
    if (
        best_family == anchor_family
        and best_identity >= identity_threshold
        and best_coverage >= coverage_threshold
    ):
        evidence_type = "seed_back_confirmed"
    return {
        "evidence_type": evidence_type,
        "identity": best_identity,
        "coverage": best_coverage,
    }


def _family_guard_terms(anchor_family: str) -> tuple[str, ...]:
    return {
        "coxL": ("coxl", "carbon monoxide", "co dehydrogenase"),
        "coxM": ("coxm", "carbon monoxide", "co dehydrogenase"),
        "coxS": ("coxs", "carbon monoxide", "co dehydrogenase"),
        "cydA": ("cyda", "appc", "cytochrome bd"),
        "cydB": ("cydb", "appb", "cytochrome bd"),
        "cydX": ("cydx",),
    }.get(anchor_family, (anchor_family.lower(),))


def _passes_family_guard(row: dict[str, object], anchor_family: str) -> bool:
    haystack = " ".join(
        str(row.get(column, ""))
        for column in ("protein_accession", "locus_tag", "gene", "product")
    ).lower()
    return any(term in haystack for term in _family_guard_terms(anchor_family))


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
    seeds_by_analyte: dict[str, dict[str, str]],
    back_confirm_identity: float,
    back_confirm_coverage: float,
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
            ],
        )
        .filter(
            pl.struct(
                ["protein_accession", "locus_tag", "gene", "product"],
            ).map_elements(
                lambda row: _passes_family_guard(row, anchor_family),
                return_dtype=pl.Boolean,
            ),
        )
        .with_columns(
            pl.struct(["sequence"])
            .map_elements(
                lambda row: _back_confirm(
                    sequence=row["sequence"],
                    analyte=analyte,
                    anchor_family=anchor_family,
                    seeds_by_analyte=seeds_by_analyte,
                    identity_threshold=back_confirm_identity,
                    coverage_threshold=back_confirm_coverage,
                ),
                return_dtype=pl.Struct(
                    {
                        "evidence_type": pl.Utf8,
                        "identity": pl.Float64,
                        "coverage": pl.Float64,
                    },
                ),
            )
            .alias("back_confirmation"),
        )
        .unnest("back_confirmation")
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


def _seed_rescue_hits(
    *,
    dataset_name: str,
    db: Path,
    analyte: str,
    anchor_family: str,
    seeds_by_analyte: dict[str, dict[str, str]],
    identity_threshold: float,
    coverage_threshold: float,
) -> pl.DataFrame:
    metadata = _all_protein_metadata(db)
    if metadata.is_empty():
        return _empty_anchor_hits()
    rows: list[dict[str, object]] = []
    seed = seeds_by_analyte.get(analyte, {}).get(anchor_family, "")
    for protein in metadata.iter_rows(named=True):
        identity, coverage = _sequence_match(str(protein["sequence"]), seed)
        if identity < identity_threshold or coverage < coverage_threshold:
            continue
        if not _passes_family_guard(protein, anchor_family):
            continue
        rows.append(
            {
                "dataset_name": dataset_name,
                "analyte": analyte,
                "anchor_family": anchor_family,
                "protein_accession": str(protein["protein_accession"]),
                "locus_tag": str(protein["locus_tag"]),
                "gene": str(protein["gene"]),
                "product": str(protein["product"]),
                "bitscore": None,
                "e_value": None,
                "identity": identity,
                "coverage": coverage,
                "evidence_type": "seed_back_confirmed",
            },
        )
    if not rows:
        return _empty_anchor_hits()
    return validate(
        pl.DataFrame(rows, schema_overrides=ANCHOR_HITS_SCHEMA),
        AnchorHitsSchema,
    )


def detect_anchors_profile(
    manifest: Path,
    *,
    config: Path | GasRegNetConfig,
    root: Path,
    profile_dir: Path,
    bitscore_threshold: float | None = None,
    e_value_threshold: float = 1e-20,
    back_confirm_identity: float = 0.25,
    back_confirm_coverage: float = 0.5,
) -> pl.DataFrame:
    """Run HMM profile anchor detection against RefSeq catalog protein FASTAs."""

    catalogs = read_refseq_catalog_manifest(manifest, root=root)
    loaded_config = load_config(config) if isinstance(config, Path) else config
    seeds_by_analyte = _seed_sequences(loaded_config)
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
                    seeds_by_analyte=seeds_by_analyte,
                    back_confirm_identity=back_confirm_identity,
                    back_confirm_coverage=back_confirm_coverage,
                )
                if not family_hits.is_empty():
                    frames.append(family_hits)
                rescue_hits = _seed_rescue_hits(
                    dataset_name=dataset_name,
                    db=db,
                    analyte=analyte.analyte,
                    anchor_family=family.name,
                    seeds_by_analyte=seeds_by_analyte,
                    identity_threshold=0.95,
                    coverage_threshold=0.95,
                )
                if not rescue_hits.is_empty():
                    frames.append(rescue_hits)
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
