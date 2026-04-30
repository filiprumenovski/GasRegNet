"""Archetype clustering.

Architecture distance is deterministic and proximity-weighted. Each architecture
is parsed into position-indexed tokens. A mismatch at position ``i`` contributes
``1 / (1 + abs(i))`` to the numerator and denominator; nearby genes therefore
matter more than distant genes. Missing positions are treated as mismatches.
"""

from __future__ import annotations

from collections import Counter
from typing import Any, cast

import polars as pl

from gasregnet.schemas import ArchetypesSchema, validate

ARCHETYPE_SCHEMA: dict[str, Any] = {
    "archetype_id": pl.Utf8,
    "analyte": pl.Utf8,
    "cluster_id": pl.Int32,
    "architecture_string": pl.Utf8,
    "n_loci": pl.Int32,
    "n_taxa": pl.Int32,
    "representative_locus_id": pl.Utf8,
    "dominant_regulator_class": pl.Utf8,
    "dominant_anchor_structure": pl.Utf8,
    "mean_locus_score": pl.Float64,
    "mean_candidate_score": pl.Float64,
}


def _candidate_token(candidate: dict[str, object]) -> str:
    relative_index = int(cast(int, candidate["relative_index"]))
    regulator_class = str(candidate["regulator_class"])
    sensory_domains = "-".join(cast(list[str], candidate["sensory_domains"])) or "none"
    return f"[{relative_index}:{regulator_class}:{sensory_domains}]"


def architecture_string(
    locus: dict[str, object],
    candidates: list[dict[str, object]],
) -> str:
    """Encode a centered architecture string for one locus."""

    candidate_tokens = sorted(
        (_candidate_token(candidate) for candidate in candidates),
        key=lambda token: int(token.split(":", 1)[0].removeprefix("[")),
    )
    anchor = f"[0:{locus['anchor_family']}]"
    return "".join([*candidate_tokens, anchor])


def _parse_architecture(architecture: str) -> dict[int, str]:
    parsed: dict[int, str] = {}
    for raw_token in architecture.split("["):
        if not raw_token:
            continue
        token = raw_token.rstrip("]")
        position_text, label = token.split(":", 1)
        parsed[int(position_text)] = label
    return parsed


def architecture_distance(left: str, right: str) -> float:
    """Return proximity-weighted normalized architecture distance."""

    left_tokens = _parse_architecture(left)
    right_tokens = _parse_architecture(right)
    positions = sorted(set(left_tokens) | set(right_tokens))
    if not positions:
        return 0.0

    numerator = 0.0
    denominator = 0.0
    for position in positions:
        weight = 1.0 / (1.0 + abs(position))
        denominator += weight
        if left_tokens.get(position) != right_tokens.get(position):
            numerator += weight
    return numerator / denominator


def _dominant(values: list[str]) -> str:
    if not values:
        return "none"
    counts = Counter(values)
    return sorted(counts.items(), key=lambda item: (-item[1], item[0]))[0][0]


def _cluster_architectures(
    locus_architectures: list[dict[str, object]],
    distance_threshold: float,
) -> list[list[dict[str, object]]]:
    clusters: list[list[dict[str, object]]] = []
    representatives: list[str] = []
    for locus_architecture in sorted(
        locus_architectures,
        key=lambda row: (str(row["architecture_string"]), str(row["locus_id"])),
    ):
        architecture = str(locus_architecture["architecture_string"])
        assigned = False
        for index, representative in enumerate(representatives):
            distance = architecture_distance(architecture, representative)
            if distance <= distance_threshold:
                clusters[index].append(locus_architecture)
                assigned = True
                break
        if not assigned:
            representatives.append(architecture)
            clusters.append([locus_architecture])
    return clusters


def cluster_archetypes(
    loci: pl.DataFrame,
    candidates: pl.DataFrame,
    *,
    distance_threshold: float = 0.0,
) -> pl.DataFrame:
    """Cluster locus architectures into deterministic archetype summaries."""

    candidates_by_locus: dict[str, list[dict[str, object]]] = {}
    for candidate in candidates.iter_rows(named=True):
        candidates_by_locus.setdefault(str(candidate["locus_id"]), []).append(candidate)

    locus_architectures: list[dict[str, object]] = []
    for locus in loci.iter_rows(named=True):
        locus_id = str(locus["locus_id"])
        locus_candidates = candidates_by_locus.get(locus_id, [])
        locus_architectures.append(
            {
                "locus_id": locus_id,
                "analyte": locus["analyte"],
                "cluster_id": locus["cluster_id"],
                "taxon_id": locus["taxon_id"],
                "anchor_family": locus["anchor_family"],
                "locus_score": locus["locus_score"],
                "candidate_score": max(
                    (
                        float(cast(float, candidate["candidate_score"]))
                        for candidate in locus_candidates
                    ),
                    default=0.0,
                ),
                "regulator_classes": [
                    str(candidate["regulator_class"]) for candidate in locus_candidates
                ],
                "architecture_string": architecture_string(locus, locus_candidates),
            },
        )

    rows: list[dict[str, object]] = []
    for index, cluster in enumerate(
        _cluster_architectures(locus_architectures, distance_threshold),
        start=1,
    ):
        representative = sorted(
            cluster,
            key=lambda row: (
                -float(cast(float, row["locus_score"])),
                str(row["locus_id"]),
            ),
        )[0]
        regulator_classes = [
            regulator_class
            for row in cluster
            for regulator_class in cast(list[str], row["regulator_classes"])
        ]
        rows.append(
            {
                "archetype_id": f"arch_{index:04d}",
                "analyte": representative["analyte"],
                "cluster_id": representative["cluster_id"],
                "architecture_string": representative["architecture_string"],
                "n_loci": len(cluster),
                "n_taxa": len({int(cast(int, row["taxon_id"])) for row in cluster}),
                "representative_locus_id": representative["locus_id"],
                "dominant_regulator_class": _dominant(regulator_classes),
                "dominant_anchor_structure": representative["anchor_family"],
                "mean_locus_score": sum(
                    float(cast(float, row["locus_score"])) for row in cluster
                )
                / len(cluster),
                "mean_candidate_score": sum(
                    float(cast(float, row["candidate_score"])) for row in cluster
                )
                / len(cluster),
            },
        )

    if not rows:
        return validate(pl.DataFrame(schema=ARCHETYPE_SCHEMA), ArchetypesSchema)
    archetypes = pl.DataFrame(rows, schema_overrides=ARCHETYPE_SCHEMA)
    return validate(archetypes, ArchetypesSchema)
