"""Ecology annotation helpers."""

from __future__ import annotations

from pathlib import Path

import polars as pl

from gasregnet.config import AnalyteConfig
from gasregnet.schemas import LociSchema, validate


def score_taxonomic_context(
    loci: pl.DataFrame,
    known_organisms: pl.DataFrame,
    *,
    matched_score: float = 1.0,
) -> pl.DataFrame:
    """Score loci whose organism or taxon appears in a known-organism table."""

    required = {"organism", "taxon_id"}
    missing = required - set(known_organisms.columns)
    if missing:
        missing_text = ", ".join(sorted(missing))
        msg = f"known-organism table is missing column(s): {missing_text}"
        raise ValueError(msg)

    known_organism_names = known_organisms["organism"].drop_nulls().cast(pl.Utf8)
    known_taxa = known_organisms["taxon_id"].drop_nulls().cast(pl.Int64)
    matched = pl.col("organism").cast(pl.Utf8).is_in(known_organism_names) | pl.col(
        "taxon_id",
    ).cast(pl.Int64).is_in(known_taxa)

    scored = loci.with_columns(
        pl.when(matched)
        .then(pl.lit(matched_score))
        .otherwise(pl.lit(0.0))
        .alias("taxonomic_context_score"),
    )
    return validate(scored, LociSchema)


def score_taxonomic_context_by_analyte(
    loci: pl.DataFrame,
    analytes: list[AnalyteConfig],
    *,
    root: Path = Path("."),
    matched_score: float = 1.0,
) -> pl.DataFrame:
    """Score loci against each analyte's configured known-organism table."""

    if loci.is_empty():
        return validate(loci, LociSchema)

    scored_parts: list[pl.DataFrame] = []
    analyte_names = [analyte.analyte for analyte in analytes]
    for analyte in analytes:
        subset = loci.filter(pl.col("analyte") == analyte.analyte)
        if subset.is_empty():
            continue
        known_path = analyte.known_organisms_table
        if not known_path.is_absolute():
            known_path = root / known_path
        known = pl.read_csv(known_path) if known_path.exists() else pl.DataFrame()
        if known.is_empty():
            scored_parts.append(
                subset.with_columns(pl.lit(0.0).alias("taxonomic_context_score")),
            )
        else:
            scored_parts.append(
                score_taxonomic_context(
                    subset,
                    known,
                    matched_score=matched_score,
                ),
            )

    missing = loci.filter(~pl.col("analyte").is_in(analyte_names))
    if not missing.is_empty():
        scored_parts.append(
            missing.with_columns(pl.lit(0.0).alias("taxonomic_context_score")),
        )
    if not scored_parts:
        return validate(loci, LociSchema)
    return validate(pl.concat(scored_parts, how="vertical"), LociSchema)
