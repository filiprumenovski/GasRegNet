"""Chemistry-partition claim helpers."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, cast

import polars as pl
from scipy.stats import chi2_contingency  # type: ignore[import-untyped]


def family_chemistry_table(candidates: pl.DataFrame) -> pl.DataFrame:
    """Count candidate regulator classes by analyte and sensory chemistry."""

    rows: list[dict[str, object]] = []
    for row in candidates.iter_rows(named=True):
        domains = cast(list[Any], row["sensory_domains"])
        if not domains:
            rows.append(
                {
                    "analyte": row["analyte"],
                    "regulator_class": row["regulator_class"],
                    "sensory_domain": "none",
                    "n_candidates": 1,
                },
            )
            continue
        for domain in domains:
            rows.append(
                {
                    "analyte": row["analyte"],
                    "regulator_class": row["regulator_class"],
                    "sensory_domain": str(domain),
                    "n_candidates": 1,
                },
            )

    if not rows:
        return pl.DataFrame(
            schema={
                "analyte": pl.Utf8,
                "regulator_class": pl.Utf8,
                "sensory_domain": pl.Utf8,
                "n_candidates": pl.Int64,
            },
        )
    return (
        pl.DataFrame(rows)
        .group_by(["analyte", "regulator_class", "sensory_domain"])
        .agg(pl.col("n_candidates").sum())
        .sort(["analyte", "regulator_class", "sensory_domain"])
    )


def chemistry_partition_outcome(
    candidates: pl.DataFrame,
    *,
    alpha: float = 0.05,
) -> dict[str, object]:
    """Test whether sensory-domain distributions differ between analytes."""

    table = family_chemistry_table(candidates)
    analytes = sorted(table["analyte"].unique().to_list()) if table.height else []
    domains = sorted(table["sensory_domain"].unique().to_list()) if table.height else []

    if len(analytes) < 2 or len(domains) < 2:
        return {
            "partition_holds": False,
            "chi_square": 0.0,
            "p_value": 1.0,
            "q_value": 1.0,
            "alpha": alpha,
            "reason": "insufficient analyte or sensory-domain diversity",
            "table": table.to_dicts(),
        }

    matrix: list[list[int]] = []
    for analyte in analytes:
        row_counts: list[int] = []
        for domain in domains:
            count = table.filter(
                (pl.col("analyte") == analyte) & (pl.col("sensory_domain") == domain),
            )["n_candidates"].sum()
            row_counts.append(int(count))
        matrix.append(row_counts)

    chi_square, p_value, _, _ = chi2_contingency(matrix)
    q_value = float(p_value)
    return {
        "partition_holds": q_value < alpha,
        "chi_square": float(chi_square),
        "p_value": float(p_value),
        "q_value": q_value,
        "alpha": alpha,
        "reason": "chi-square test on analyte-by-sensory-domain counts",
        "analytes": analytes,
        "sensory_domains": domains,
        "table": table.to_dicts(),
    }


def write_partition_outcome(outcome: dict[str, object], out_path: Path) -> Path:
    """Write partition outcome JSON."""

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as handle:
        json.dump(outcome, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return out_path
