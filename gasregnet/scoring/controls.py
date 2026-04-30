"""Matched-control sampling."""

from __future__ import annotations

import random
from typing import Any

import polars as pl

from gasregnet.schemas import LociSchema, validate

LOCI_SCHEMA_OVERRIDES: dict[str, Any] = {
    "cluster_id": pl.Int32,
    "window_size": pl.Int32,
    "marker_genes_present": pl.List(pl.Utf8),
    "accessory_genes_present": pl.List(pl.Utf8),
    "created_at": pl.Datetime("us"),
}


def _rows_by_match_key(
    rows: list[dict[str, object]],
) -> dict[tuple[str, str], list[int]]:
    grouped: dict[tuple[str, str], list[int]] = {}
    for index, row in enumerate(rows):
        key = (str(row["organism"]), str(row["contig_id"]))
        grouped.setdefault(key, []).append(index)
    return grouped


def sample_matched_controls(
    case_loci: pl.DataFrame,
    candidate_pool: pl.DataFrame,
    *,
    ratio: tuple[int, int],
    seed: int,
) -> pl.DataFrame:
    """Sample control loci matched on organism and contig."""

    case_rows = case_loci.to_dicts()
    pool_rows = [
        row
        for row in candidate_pool.to_dicts()
        if row["locus_confidence"] != "high"
        and row["locus_id"] not in set(case_loci["locus_id"].to_list())
    ]
    grouped_pool = _rows_by_match_key(pool_rows)

    numerator, denominator = ratio
    if numerator <= 0 or denominator < 0:
        msg = f"ratio must be positive case:control values, got {ratio}"
        raise ValueError(msg)

    rng = random.Random(seed)
    selected: list[dict[str, object]] = []
    used_indices: set[int] = set()
    per_case = denominator / numerator

    for case_row in sorted(case_rows, key=lambda row: str(row["locus_id"])):
        key = (str(case_row["organism"]), str(case_row["contig_id"]))
        candidate_indices = [
            index for index in grouped_pool.get(key, []) if index not in used_indices
        ]
        n_controls = min(int(per_case), len(candidate_indices))
        sampled_indices = rng.sample(candidate_indices, n_controls)
        for index in sampled_indices:
            used_indices.add(index)
            control_row = dict(pool_rows[index])
            control_row["locus_confidence"] = "control"
            selected.append(control_row)

    if not selected:
        return validate(candidate_pool.head(0), LociSchema)

    controls = pl.DataFrame(selected, schema_overrides=LOCI_SCHEMA_OVERRIDES)
    return validate(controls, LociSchema)
