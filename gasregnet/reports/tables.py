"""Publication table builders."""

from __future__ import annotations

from pathlib import Path

import polars as pl

TABLE_SPECS = {
    "T1_benchmark_recovery": [
        "benchmark_id",
        "analyte",
        "protein_name",
        "organism",
        "hit",
        "rank",
        "candidate_score",
    ],
    "T2_top30_co_candidates": [
        "candidate_id",
        "organism",
        "regulator_class",
        "sensory_domains",
        "candidate_score",
        "candidate_score_q",
    ],
    "T3_top30_hcn_candidates": [
        "candidate_id",
        "organism",
        "regulator_class",
        "sensory_domains",
        "candidate_score",
        "candidate_score_q",
    ],
    "T4_regulator_family_enrichment": [
        "analyte",
        "feature_type",
        "feature_name",
        "odds_ratio",
        "p_value",
        "q_value",
        "interpretation",
    ],
    "T5_archetype_catalog": [
        "archetype_id",
        "analyte",
        "architecture_string",
        "n_loci",
        "n_taxa",
        "dominant_regulator_class",
    ],
    "T6_tool_feature_comparison": [
        "tool",
        "analyte_general",
        "decomposable_scoring",
        "matched_control_enrichment",
        "archetype_clustering",
        "structural_prioritization",
        "fdr_controlled_candidates",
        "reproducible_workflow",
    ],
}


def _select_existing(frame: pl.DataFrame, columns: list[str]) -> pl.DataFrame:
    existing = [column for column in columns if column in frame.columns]
    return frame.select(existing) if existing else frame


def _stringify_nested(frame: pl.DataFrame) -> pl.DataFrame:
    expressions = []
    for column, dtype in zip(frame.columns, frame.dtypes, strict=True):
        if dtype.base_type() is pl.List:
            expressions.append(pl.col(column).list.join(";").alias(column))
        else:
            expressions.append(pl.col(column))
    return frame.select(expressions)


def _write_table(frame: pl.DataFrame, stem: str, out_dir: Path) -> tuple[Path, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / f"{stem}.csv"
    markdown_path = out_dir / f"{stem}.md"
    frame = _stringify_nested(frame)
    frame.write_csv(csv_path)
    markdown_path.write_text(
        frame.to_pandas().to_markdown(index=False) + "\n",
        encoding="utf-8",
    )
    return csv_path, markdown_path


def tool_feature_comparison() -> pl.DataFrame:
    """Return the manuscript tool-comparison matrix."""

    return pl.DataFrame(
        {
            "tool": ["EFI-GNT", "RODEO", "antiSMASH", "FlaGs/webFlaGs", "GasRegNet"],
            "analyte_general": [False, False, False, False, True],
            "decomposable_scoring": [False, False, False, False, True],
            "matched_control_enrichment": [False, False, False, False, True],
            "archetype_clustering": [False, True, True, False, True],
            "structural_prioritization": [False, False, False, False, True],
            "fdr_controlled_candidates": [False, False, False, False, True],
            "reproducible_workflow": [False, False, True, False, True],
        },
    )


def write_publication_tables(
    *,
    benchmark_recovery: pl.DataFrame,
    candidates: pl.DataFrame,
    enrichment: pl.DataFrame,
    archetypes: pl.DataFrame,
    out_dir: Path,
) -> dict[str, tuple[Path, Path]]:
    """Write T1-T6 as CSV and Markdown files."""

    tables = {
        "T1_benchmark_recovery": benchmark_recovery,
        "T2_top30_co_candidates": candidates.filter(pl.col("analyte") == "CO")
        .sort("candidate_score", descending=True)
        .head(30),
        "T3_top30_hcn_candidates": candidates.filter(pl.col("analyte") == "CN")
        .sort("candidate_score", descending=True)
        .head(30),
        "T4_regulator_family_enrichment": enrichment,
        "T5_archetype_catalog": archetypes,
        "T6_tool_feature_comparison": tool_feature_comparison(),
    }
    outputs: dict[str, tuple[Path, Path]] = {}
    for stem, frame in tables.items():
        outputs[stem] = _write_table(
            _select_existing(frame, TABLE_SPECS[stem]),
            stem,
            out_dir,
        )
    return outputs
