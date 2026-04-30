"""Matched-control enrichment statistics."""

from __future__ import annotations

from typing import Any, cast

import polars as pl
from scipy.stats import fisher_exact  # type: ignore[import-untyped]
from statsmodels.stats.multitest import multipletests  # type: ignore[import-untyped]

from gasregnet.schemas import EnrichmentResultsSchema, validate

ENRICHMENT_SCHEMA_OVERRIDES: dict[str, Any] = {
    "n_case_with_feature": pl.Int64,
    "n_case_without_feature": pl.Int64,
    "n_control_with_feature": pl.Int64,
    "n_control_without_feature": pl.Int64,
    "odds_ratio": pl.Float64,
    "p_value": pl.Float64,
    "q_value": pl.Float64,
}
ENRICHMENT_EMPTY_SCHEMA: dict[str, Any] = {
    "analyte": pl.Utf8,
    "feature_type": pl.Utf8,
    "feature_name": pl.Utf8,
    "case_definition": pl.Utf8,
    "control_definition": pl.Utf8,
    "n_case_with_feature": pl.Int64,
    "n_case_without_feature": pl.Int64,
    "n_control_with_feature": pl.Int64,
    "n_control_without_feature": pl.Int64,
    "odds_ratio": pl.Float64,
    "p_value": pl.Float64,
    "q_value": pl.Float64,
    "interpretation": pl.Utf8,
}


def _feature_loci(genes: pl.DataFrame, feature_type: str) -> dict[str, set[str]]:
    feature_loci: dict[str, set[str]] = {}
    if feature_type == "regulator_class":
        rows = genes.filter(
            (pl.col("is_regulator_candidate")) & (pl.col("regulator_class") != "none"),
        ).select(["locus_id", "regulator_class"])
        for row in rows.iter_rows(named=True):
            feature_loci.setdefault(str(row["regulator_class"]), set()).add(
                str(row["locus_id"]),
            )
        return feature_loci

    if feature_type == "sensory_domain":
        rows = (
            genes.select(["locus_id", "sensory_domains"])
            .explode("sensory_domains")
            .drop_nulls("sensory_domains")
            .filter(pl.col("sensory_domains") != "")
        )
        for row in rows.iter_rows(named=True):
            feature_loci.setdefault(str(row["sensory_domains"]), set()).add(
                str(row["locus_id"]),
            )
        return feature_loci

    msg = f"unsupported feature type: {feature_type}"
    raise ValueError(msg)


def _total_loci(genes: pl.DataFrame) -> set[str]:
    return {str(locus_id) for locus_id in genes["locus_id"].unique().to_list()}


def _rows_for_feature_type(
    case_genes: pl.DataFrame,
    control_genes: pl.DataFrame,
    *,
    analyte: str,
    feature_type: str,
    case_definition: str,
    control_definition: str,
) -> list[dict[str, object]]:
    case_loci = _total_loci(case_genes)
    control_loci = _total_loci(control_genes)
    case_features = _feature_loci(case_genes, feature_type)
    control_features = _feature_loci(control_genes, feature_type)
    feature_names = sorted(set(case_features) | set(control_features))

    rows: list[dict[str, object]] = []
    for feature_name in feature_names:
        case_with = len(case_features.get(feature_name, set()))
        control_with = len(control_features.get(feature_name, set()))
        case_without = len(case_loci) - case_with
        control_without = len(control_loci) - control_with
        odds_ratio, p_value = fisher_exact(
            [[case_with, case_without], [control_with, control_without]],
            alternative="two-sided",
        )
        rows.append(
            {
                "analyte": analyte,
                "feature_type": feature_type,
                "feature_name": feature_name,
                "case_definition": case_definition,
                "control_definition": control_definition,
                "n_case_with_feature": case_with,
                "n_case_without_feature": case_without,
                "n_control_with_feature": control_with,
                "n_control_without_feature": control_without,
                "odds_ratio": float(odds_ratio),
                "p_value": float(p_value),
                "q_value": 1.0,
                "interpretation": "not_tested",
            },
        )
    return rows


def run_enrichment(
    case_genes: pl.DataFrame,
    control_genes: pl.DataFrame,
    *,
    analyte: str,
    case_definition: str,
    control_definition: str,
    alpha: float = 0.05,
) -> pl.DataFrame:
    """Run Fisher exact enrichment with BH correction per analyte."""

    rows: list[dict[str, object]] = []
    for feature_type in ("regulator_class", "sensory_domain"):
        rows.extend(
            _rows_for_feature_type(
                case_genes,
                control_genes,
                analyte=analyte,
                feature_type=feature_type,
                case_definition=case_definition,
                control_definition=control_definition,
            ),
    )

    if not rows:
        return validate(
            pl.DataFrame(schema=ENRICHMENT_EMPTY_SCHEMA),
            EnrichmentResultsSchema,
        )

    p_values = [float(cast(float, row["p_value"])) for row in rows]
    _, q_values, _, _ = multipletests(p_values, alpha=alpha, method="fdr_bh")
    for row, q_value in zip(rows, q_values, strict=True):
        row["q_value"] = float(q_value)
        row["interpretation"] = "enriched" if q_value <= alpha else "not_enriched"

    enrichment = pl.DataFrame(rows, schema_overrides=ENRICHMENT_SCHEMA_OVERRIDES)
    return validate(enrichment, EnrichmentResultsSchema)
