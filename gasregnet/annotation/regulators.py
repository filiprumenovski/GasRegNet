"""Regulator classification helpers."""

from __future__ import annotations

from typing import Any

import polars as pl

from gasregnet.config import RegulatorFamilyEntry
from gasregnet.schemas import GenesSchema, validate

ALLOWED_NON_REGULATOR_CLASSES = {"anchor", "metabolic", "transporter", "unknown"}
GENES_SCHEMA_OVERRIDES: dict[str, Any] = {
    "relative_index": pl.Int32,
    "relative_start": pl.Int64,
    "relative_stop": pl.Int64,
    "pfam_ids": pl.List(pl.Utf8),
    "pfam_descriptions": pl.List(pl.Utf8),
    "interpro_ids": pl.List(pl.Utf8),
    "interpro_descriptions": pl.List(pl.Utf8),
    "sensory_domains": pl.List(pl.Utf8),
    "is_anchor": pl.Boolean,
    "is_regulator_candidate": pl.Boolean,
}


def _has_any(pfam_ids: set[str], required: set[str]) -> bool:
    return bool(pfam_ids & required)


def _has_all(pfam_ids: set[str], required: set[str]) -> bool:
    return required.issubset(pfam_ids)


def _classify_regulator(
    pfam_ids: set[str],
    regulator_families: list[RegulatorFamilyEntry],
) -> str:
    scores: dict[str, tuple[float, int]] = {}
    if "PF00072" in pfam_ids:
        scores["two_component_rr"] = (1.0, 10)
    if _has_any(pfam_ids, {"PF00512", "PF02518"}):
        scores["two_component_hk"] = (1.0, 9)
    if _has_any(pfam_ids, {"PF04542", "PF04545"}):
        scores["sigma"] = (1.0, 8)
    for entry in regulator_families:
        if _has_all(pfam_ids, set(entry.pfam_required)):
            score = 2.0 + 0.5 * len(set(entry.pfam_optional) & pfam_ids)
            priority = 100 - regulator_families.index(entry)
            current = scores.get(entry.regulator_class)
            if current is None or (score, priority) > current:
                scores[entry.regulator_class] = (score, priority)
    if not scores:
        return "none"
    return max(scores.items(), key=lambda item: (item[1][0], item[1][1]))[0]
    return "none"


def classify_regulators(
    genes: pl.DataFrame,
    regulator_families: list[RegulatorFamilyEntry],
) -> pl.DataFrame:
    """Classify regulatory genes from Pfam rules in configuration."""

    rows: list[dict[str, object]] = []
    for row in genes.iter_rows(named=True):
        updated = dict(row)
        pfam_ids = {str(pfam_id) for pfam_id in row["pfam_ids"]}
        regulator_class = _classify_regulator(pfam_ids, regulator_families)
        is_candidate = regulator_class != "none"
        updated["regulator_class"] = regulator_class
        updated["is_regulator_candidate"] = is_candidate

        if bool(row["is_anchor"]):
            updated["functional_class"] = "anchor"
        elif is_candidate:
            updated["functional_class"] = "regulator"
        elif row["functional_class"] in ALLOWED_NON_REGULATOR_CLASSES:
            updated["functional_class"] = row["functional_class"]
        else:
            updated["functional_class"] = "unknown"
        rows.append(updated)

    classified = pl.DataFrame(
        rows,
        schema_overrides=GENES_SCHEMA_OVERRIDES,
    )
    return validate(classified, GenesSchema)
