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

    known_parts: list[pl.DataFrame] = []
    for analyte in analytes:
        known_path = analyte.known_organisms_table
        if not known_path.is_absolute():
            known_path = root / known_path
        known = pl.read_csv(known_path) if known_path.exists() else pl.DataFrame()
        if known.is_empty():
            continue
        required = {"organism", "taxon_id"}
        missing_columns = required - set(known.columns)
        if missing_columns:
            missing_text = ", ".join(sorted(missing_columns))
            msg = f"known-organism table is missing column(s): {missing_text}"
            raise ValueError(msg)
        known_parts.append(
            known.select(
                pl.lit(analyte.analyte).alias("analyte"),
                pl.col("organism").cast(pl.Utf8),
                pl.col("taxon_id").cast(pl.Int64),
            )
            .unique(maintain_order=True),
        )

    if not known_parts:
        return validate(
            loci.with_columns(pl.lit(0.0).alias("taxonomic_context_score")),
            LociSchema,
        )

    known = pl.concat(known_parts, how="vertical")
    organism_matches = known.select("analyte", "organism").unique().with_columns(
        pl.lit(True).alias("__organism_match"),
    )
    taxon_matches = known.select("analyte", "taxon_id").unique().with_columns(
        pl.lit(True).alias("__taxon_match"),
    )
    analyte_order = pl.DataFrame(
        {
            "analyte": [analyte.analyte for analyte in analytes],
            "__analyte_order": list(range(len(analytes))),
        },
    )
    scored = (
        loci.with_row_index("__row_index")
        .join(organism_matches, on=["analyte", "organism"], how="left")
        .join(taxon_matches, on=["analyte", "taxon_id"], how="left")
        .join(analyte_order, on="analyte", how="left")
        .with_columns(
            pl.when(
                pl.col("__organism_match").fill_null(False)
                | pl.col("__taxon_match").fill_null(False),
            )
            .then(pl.lit(matched_score))
            .otherwise(pl.lit(0.0))
            .alias("taxonomic_context_score"),
            pl.col("__analyte_order").fill_null(len(analytes)),
        )
        .sort(["__analyte_order", "__row_index"])
        .drop(
            [
                "__row_index",
                "__organism_match",
                "__taxon_match",
                "__analyte_order",
            ],
        )
    )
    return validate(scored, LociSchema)
