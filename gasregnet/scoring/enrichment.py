"""Matched-control enrichment statistics."""

from __future__ import annotations

import math
from typing import Any, cast

import polars as pl
from scipy.stats import fisher_exact  # type: ignore[import-untyped]
from statsmodels.stats.contingency_tables import (  # type: ignore[import-untyped]
    StratifiedTable,
)
from statsmodels.stats.multitest import multipletests  # type: ignore[import-untyped]

from gasregnet.schemas import (
    EnrichmentResultsSchema,
    EnrichmentRobustnessSchema,
    validate,
)

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
ROBUSTNESS_SCHEMA_OVERRIDES: dict[str, Any] = {
    "q_value": pl.Float64,
    "p_value": pl.Float64,
    "odds_ratio": pl.Float64,
}
ROBUSTNESS_EMPTY_SCHEMA: dict[str, Any] = {
    "analyte": pl.Utf8,
    "feature_type": pl.Utf8,
    "feature_name": pl.Utf8,
    "deduplication_policy": pl.Utf8,
    "stratum_column": pl.Utf8,
    "test": pl.Utf8,
    "q_value": pl.Float64,
    "p_value": pl.Float64,
    "odds_ratio": pl.Float64,
}


def _taxonomy_value(row: dict[str, object], column: str) -> str:
    value = row.get(column)
    if isinstance(value, str) and value:
        return value
    organism = row.get("organism")
    if column == "genus" and isinstance(organism, str) and organism:
        return organism.split()[0]
    if column == "family":
        taxon_id = row.get("taxon_id")
        if taxon_id is not None:
            return str(taxon_id)
    return "unknown"


def _feature_values(row: dict[str, object], feature_type: str) -> set[str]:
    if feature_type == "regulator_class":
        if (
            bool(row.get("is_regulator_candidate"))
            and row.get("regulator_class") != "none"
        ):
            return {str(row["regulator_class"])}
        return set()
    if feature_type == "sensory_domain":
        domains = row.get("sensory_domains")
        if isinstance(domains, list):
            return {str(domain) for domain in domains if domain}
        return set()
    msg = f"unsupported feature type: {feature_type}"
    raise ValueError(msg)


def _locus_records(
    genes: pl.DataFrame,
    *,
    policy: str = "none",
    stratum_column: str = "genus",
) -> dict[str, dict[str, object]]:
    records: dict[str, dict[str, object]] = {}
    for row in genes.iter_rows(named=True):
        locus_id = str(row["locus_id"])
        record = records.setdefault(
            locus_id,
            {
                "locus_id": locus_id,
                "genus": _taxonomy_value(row, "genus"),
                "family": _taxonomy_value(row, "family"),
                "stratum": _taxonomy_value(row, stratum_column),
                "locus_score": float(cast(float, row.get("locus_score", 0.0) or 0.0)),
                "regulator_class": set(),
                "sensory_domain": set(),
            },
        )
        cast(set[str], record["regulator_class"]).update(
            _feature_values(row, "regulator_class"),
        )
        cast(set[str], record["sensory_domain"]).update(
            _feature_values(row, "sensory_domain"),
        )

    if policy == "none":
        return records
    if policy not in {"one_per_genus", "one_per_family"}:
        msg = f"unsupported deduplication policy: {policy}"
        raise ValueError(msg)

    group_column = "genus" if policy == "one_per_genus" else "family"
    selected: dict[str, dict[str, object]] = {}
    for record in records.values():
        group = str(record[group_column])
        current = selected.get(group)
        if current is None or float(cast(float, record["locus_score"])) > float(
            cast(float, current["locus_score"]),
        ):
            selected[group] = record
    return {str(record["locus_id"]): record for record in selected.values()}


def _feature_loci_from_records(
    records: dict[str, dict[str, object]],
    feature_type: str,
) -> dict[str, set[str]]:
    feature_loci: dict[str, set[str]] = {}
    for locus_id, record in records.items():
        for feature_name in cast(set[str], record[feature_type]):
            feature_loci.setdefault(feature_name, set()).add(locus_id)
    return feature_loci


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
        if math.isnan(odds_ratio):
            odds_ratio = 0.0
        if math.isnan(p_value):
            p_value = 1.0
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


def _cmh_for_feature(
    case_records: dict[str, dict[str, object]],
    control_records: dict[str, dict[str, object]],
    *,
    feature_type: str,
    feature_name: str,
) -> tuple[int, int, int, int, float, float]:
    case_loci = set(case_records)
    control_loci = set(control_records)
    case_features = _feature_loci_from_records(case_records, feature_type)
    control_features = _feature_loci_from_records(control_records, feature_type)
    case_with_loci = case_features.get(feature_name, set())
    control_with_loci = control_features.get(feature_name, set())
    case_without_loci = case_loci - case_with_loci
    control_without_loci = control_loci - control_with_loci

    strata = sorted(
        {str(record["stratum"]) for record in case_records.values()}
        | {str(record["stratum"]) for record in control_records.values()},
    )
    tables: list[list[list[int]]] = []
    for stratum in strata:
        case_in_stratum = {
            locus_id
            for locus_id, record in case_records.items()
            if str(record["stratum"]) == stratum
        }
        control_in_stratum = {
            locus_id
            for locus_id, record in control_records.items()
            if str(record["stratum"]) == stratum
        }
        table = [
            [
                len(case_with_loci & case_in_stratum),
                len(case_without_loci & case_in_stratum),
            ],
            [
                len(control_with_loci & control_in_stratum),
                len(control_without_loci & control_in_stratum),
            ],
        ]
        if sum(sum(row) for row in table) > 0:
            tables.append(table)

    if not tables:
        return (0, 0, 0, 0, 0.0, 1.0)

    try:
        stratified = StratifiedTable(tables, shift_zeros=True)
        result = stratified.test_null_odds(correction=True)
        odds_ratio = float(stratified.oddsratio_pooled)
        p_value = float(result.pvalue)
    except Exception:
        odds_ratio, p_value = fisher_exact(
            [
                [len(case_with_loci), len(case_without_loci)],
                [len(control_with_loci), len(control_without_loci)],
            ],
            alternative="two-sided",
        )
    if math.isnan(odds_ratio):
        odds_ratio = 0.0
    if math.isnan(p_value):
        p_value = 1.0
    return (
        len(case_with_loci),
        len(case_without_loci),
        len(control_with_loci),
        len(control_without_loci),
        odds_ratio,
        p_value,
    )


def _stratified_rows_for_feature_type(
    case_genes: pl.DataFrame,
    control_genes: pl.DataFrame,
    *,
    analyte: str,
    feature_type: str,
    case_definition: str,
    control_definition: str,
    stratum_column: str,
    deduplication_policy: str,
) -> list[dict[str, object]]:
    case_records = _locus_records(
        case_genes,
        policy=deduplication_policy,
        stratum_column=stratum_column,
    )
    control_records = _locus_records(
        control_genes,
        policy=deduplication_policy,
        stratum_column=stratum_column,
    )
    case_features = _feature_loci_from_records(case_records, feature_type)
    control_features = _feature_loci_from_records(control_records, feature_type)
    feature_names = sorted(set(case_features) | set(control_features))

    rows: list[dict[str, object]] = []
    for feature_name in feature_names:
        (
            case_with,
            case_without,
            control_with,
            control_without,
            odds_ratio,
            p_value,
        ) = _cmh_for_feature(
            case_records,
            control_records,
            feature_type=feature_type,
            feature_name=feature_name,
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


def run_stratified_enrichment(
    case_genes: pl.DataFrame,
    control_genes: pl.DataFrame,
    *,
    analyte: str,
    case_definition: str,
    control_definition: str,
    stratum_column: str = "genus",
    test: str = "cmh",
    alpha: float = 0.05,
    deduplication_policy: str = "one_per_family",
) -> pl.DataFrame:
    """Run CMH enrichment with BH correction across tested features."""

    if test != "cmh":
        return run_enrichment(
            case_genes,
            control_genes,
            analyte=analyte,
            case_definition=case_definition,
            control_definition=control_definition,
            alpha=alpha,
        )

    rows: list[dict[str, object]] = []
    for feature_type in ("regulator_class", "sensory_domain"):
        rows.extend(
            _stratified_rows_for_feature_type(
                case_genes,
                control_genes,
                analyte=analyte,
                feature_type=feature_type,
                case_definition=case_definition,
                control_definition=control_definition,
                stratum_column=stratum_column,
                deduplication_policy=deduplication_policy,
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


def run_enrichment_robustness(
    case_genes: pl.DataFrame,
    control_genes: pl.DataFrame,
    *,
    analyte: str,
    case_definition: str,
    control_definition: str,
    stratum_column: str = "genus",
    alpha: float = 0.05,
) -> pl.DataFrame:
    """Run CMH enrichment under no, genus, and family deduplication policies."""

    rows: list[dict[str, object]] = []
    for policy in ("none", "one_per_genus", "one_per_family"):
        enrichment = run_stratified_enrichment(
            case_genes,
            control_genes,
            analyte=analyte,
            case_definition=case_definition,
            control_definition=control_definition,
            stratum_column=stratum_column,
            test="cmh",
            alpha=alpha,
            deduplication_policy=policy,
        )
        for row in enrichment.iter_rows(named=True):
            rows.append(
                {
                    "analyte": row["analyte"],
                    "feature_type": row["feature_type"],
                    "feature_name": row["feature_name"],
                    "deduplication_policy": policy,
                    "stratum_column": stratum_column,
                    "test": "cmh",
                    "q_value": row["q_value"],
                    "p_value": row["p_value"],
                    "odds_ratio": row["odds_ratio"],
                },
            )

    if not rows:
        return validate(
            pl.DataFrame(schema=ROBUSTNESS_EMPTY_SCHEMA),
            EnrichmentRobustnessSchema,
        )
    robustness = pl.DataFrame(rows, schema_overrides=ROBUSTNESS_SCHEMA_OVERRIDES)
    return validate(robustness, EnrichmentRobustnessSchema)
