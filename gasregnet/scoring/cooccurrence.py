"""Phylogenetic profile co-occurrence evidence for regulator candidates."""

from __future__ import annotations

import math
from typing import Any, cast

import polars as pl

from gasregnet.config import ScoringConfig
from gasregnet.schemas import RegulatorCandidatesSchema, validate
from gasregnet.scoring.candidates import CANDIDATE_SCHEMA


def _entropy(probability: float) -> float:
    if probability <= 0.0 or probability >= 1.0:
        return 0.0
    return -(
        probability * math.log2(probability)
        + (1.0 - probability) * math.log2(1.0 - probability)
    )


def _mutual_information(x: list[int], y: list[int]) -> float:
    n = len(x)
    if n == 0:
        return 0.0
    mi = 0.0
    for x_value in (0, 1):
        px = sum(1 for value in x if value == x_value) / n
        if px == 0.0:
            continue
        for y_value in (0, 1):
            py = sum(1 for value in y if value == y_value) / n
            pxy = (
                sum(
                    1
                    for left, right in zip(x, y, strict=True)
                    if left == x_value and right == y_value
                )
                / n
            )
            if pxy > 0.0 and py > 0.0:
                mi += pxy * math.log2(pxy / (px * py))
    normalizer = min(_entropy(sum(x) / n), _entropy(sum(y) / n))
    if normalizer == 0.0:
        return 0.0
    return max(0.0, min(1.0, mi / normalizer))


def _score_total(row: dict[str, object], scoring: ScoringConfig) -> float:
    weights = scoring.candidate_score_weights
    structural_score = row["structural_plausibility_score"]
    structural_component = (
        0.0 if structural_score is None else float(cast(float, structural_score))
    )
    return (
        float(cast(float, row["locus_score"])) * weights["locus"]
        + float(cast(float, row["regulator_domain_score"]))
        * weights["regulator_domain"]
        + float(cast(float, row["sensory_domain_score"])) * weights["sensory_domain"]
        + float(cast(float, row["proximity_score"])) * weights["proximity"]
        + float(cast(float, row["archetype_conservation_score"]))
        * weights["archetype_conservation"]
        + float(cast(float, row["enrichment_score"])) * weights["enrichment"]
        + float(cast(float, row["taxonomic_breadth_score"]))
        * weights.get("taxonomic_breadth", 0.0)
        + float(cast(float, row.get("phylogenetic_profile_score", 0.0)))
        * weights.get("phylogenetic_profile", 0.0)
        + structural_component * weights["structural_plausibility"]
    )


def assign_phylogenetic_profile_scores(
    candidates: pl.DataFrame,
    loci: pl.DataFrame,
    *,
    scoring: ScoringConfig | None = None,
) -> pl.DataFrame:
    """Attach normalized mutual-information scores across organisms.

    The profile is computed per ``(analyte, regulator_class)`` as regulator-class
    presence/absence versus anchor presence/absence across organisms. Single-genome
    runs receive score 0.0 by construction because no phylogenetic profile exists.
    """

    candidates = validate(candidates, RegulatorCandidatesSchema)
    if candidates.is_empty() or loci.is_empty():
        return candidates

    organisms = sorted(
        {
            str(value)
            for value in loci["organism"].to_list() + candidates["organism"].to_list()
            if value is not None
        },
    )
    if len(organisms) < 2:
        if "phylogenetic_profile_score" in candidates.columns:
            return candidates
        return validate(
            candidates.with_columns(pl.lit(0.0).alias("phylogenetic_profile_score")),
            RegulatorCandidatesSchema,
        )

    anchor_presence: dict[tuple[str, str], set[str]] = {}
    for locus in loci.iter_rows(named=True):
        anchor_presence.setdefault(
            (str(locus["analyte"]), str(locus["organism"])),
            set(),
        ).add(str(locus["anchor_family"]))

    regulator_presence: dict[tuple[str, str, str], bool] = {}
    for row in candidates.iter_rows(named=True):
        regulator_presence[
            (str(row["analyte"]), str(row["regulator_class"]), str(row["organism"]))
        ] = True

    score_by_key: dict[tuple[str, str], float] = {}
    for analyte, regulator_class in (
        candidates.select(
            ["analyte", "regulator_class"],
        )
        .unique()
        .iter_rows()
    ):
        anchor_vector = [
            int(bool(anchor_presence.get((str(analyte), organism))))
            for organism in organisms
        ]
        regulator_vector = [
            int(
                regulator_presence.get(
                    (str(analyte), str(regulator_class), organism),
                    False,
                ),
            )
            for organism in organisms
        ]
        score_by_key[(str(analyte), str(regulator_class))] = _mutual_information(
            anchor_vector,
            regulator_vector,
        )

    rows: list[dict[str, Any]] = []
    for row in candidates.iter_rows(named=True):
        updated = dict(row)
        updated["phylogenetic_profile_score"] = score_by_key.get(
            (str(row["analyte"]), str(row["regulator_class"])),
            0.0,
        )
        if scoring is not None:
            updated["candidate_score"] = _score_total(updated, scoring)
        rows.append(updated)

    return validate(
        pl.DataFrame(rows, schema_overrides=CANDIDATE_SCHEMA),
        RegulatorCandidatesSchema,
    )
