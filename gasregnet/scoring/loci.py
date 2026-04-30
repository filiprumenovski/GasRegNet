"""Locus scoring."""

from __future__ import annotations

from typing import Any

import polars as pl

from gasregnet.config import ScoringConfig
from gasregnet.schemas import LociSchema, validate

SCORED_LOCI_SCHEMA_OVERRIDES: dict[str, Any] = {
    "cluster_id": pl.Int32,
    "window_size": pl.Int32,
    "marker_genes_present": pl.List(pl.Utf8),
    "accessory_genes_present": pl.List(pl.Utf8),
    "created_at": pl.Datetime("us"),
}


def anchor_marker_score() -> pl.Expr:
    """Score direct gas-marker evidence from marker gene presence."""

    return (
        pl.when(pl.col("marker_genes_present").list.len() > 0)
        .then(1.0)
        .otherwise(0.0)
    )


def accessory_marker_score() -> pl.Expr:
    """Score accessory-marker support as a capped marker count."""

    return pl.min_horizontal(
        pl.col("accessory_genes_present").list.len().cast(pl.Float64),
        pl.lit(3.0),
    )


def operon_integrity_component_score() -> pl.Expr:
    """Use the existing operon-integrity estimate as a scoring component."""

    return pl.col("operon_integrity_score").cast(pl.Float64)


def homology_confidence_score() -> pl.Expr:
    """Score loci with a non-empty anchor accession as having homology support."""

    return (
        pl.when(pl.col("anchor_accession").str.len_chars() > 0)
        .then(1.0)
        .otherwise(0.0)
    )


def taxonomic_context_component_score() -> pl.Expr:
    """Use the existing taxonomic-context estimate as a scoring component."""

    return pl.col("taxonomic_context_score").cast(pl.Float64)


def neighborhood_completeness_score() -> pl.Expr:
    """Score complete windows higher than boundary-truncated windows."""

    return pl.when(pl.col("is_boundary_truncated")).then(0.0).otherwise(1.0)


def _confidence_expr(scoring: ScoringConfig) -> pl.Expr:
    thresholds = scoring.confidence_thresholds
    return (
        pl.when(pl.col("locus_score") >= thresholds.high)
        .then(pl.lit("high"))
        .when(pl.col("locus_score") >= thresholds.medium)
        .then(pl.lit("medium"))
        .when(pl.col("locus_score") >= thresholds.low)
        .then(pl.lit("low"))
        .otherwise(pl.lit("control"))
    )


def score_loci(loci: pl.DataFrame, scoring: ScoringConfig) -> pl.DataFrame:
    """Compute deterministic decomposable locus scores."""

    weights = scoring.locus_score_weights
    scored = loci.with_columns(
        anchor_marker_score().alias("anchor_marker_score"),
        accessory_marker_score().alias("accessory_marker_score"),
        operon_integrity_component_score().alias("operon_integrity_component_score"),
        homology_confidence_score().alias("homology_confidence_score"),
        taxonomic_context_component_score().alias("taxonomic_context_component_score"),
        neighborhood_completeness_score().alias("neighborhood_completeness_score"),
    ).with_columns(
        (
            pl.col("anchor_marker_score") * weights["anchor_marker"]
            + pl.col("accessory_marker_score") * weights["accessory_marker"]
            + pl.col("operon_integrity_component_score")
            * weights["operon_integrity"]
            + pl.col("homology_confidence_score") * weights["homology_confidence"]
            + pl.col("taxonomic_context_component_score")
            * weights["taxonomic_context"]
            + pl.col("neighborhood_completeness_score")
            * weights["neighborhood_completeness"]
        ).alias("locus_score"),
    )
    scored = scored.with_columns(_confidence_expr(scoring).alias("locus_confidence"))

    base_columns = list(LociSchema.columns.keys())
    validate(scored.select(base_columns), LociSchema)
    return pl.DataFrame(scored, schema_overrides=SCORED_LOCI_SCHEMA_OVERRIDES)
