"""Publication table builders."""

from __future__ import annotations

from pathlib import Path

import polars as pl

from gasregnet.benchmark import summarize_benchmark_recovery

TABLE_SPECS = {
    "T0_benchmark_validation_summary": [
        "metric",
        "value",
        "n",
        "notes",
    ],
    "T1_benchmark_recovery": [
        "benchmark_id",
        "analyte",
        "sensing_evidence_class",
        "is_negative_control",
        "verified_pmid",
        "protein_name",
        "organism",
        "hit",
        "rank",
        "regulation_logit_score",
        "candidate_score",
    ],
    "T2_top30_co_candidates": [
        "candidate_id",
        "organism",
        "regulator_class",
        "sensory_domains",
        "regulation_logit_score",
        "score_band_low",
        "score_band_high",
        "phylogenetic_profile_score",
        "candidate_score",
        "candidate_score_q",
    ],
    "T3_top30_cyd_control_candidates": [
        "candidate_id",
        "organism",
        "regulator_class",
        "sensory_domains",
        "regulation_logit_score",
        "score_band_low",
        "score_band_high",
        "phylogenetic_profile_score",
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
        "architecture_scope",
        "architecture_string",
        "n_loci",
        "n_taxa",
        "dominant_regulator_class",
    ],
    "T6_tool_feature_comparison": [
        "tool",
        "analyte_general",
        "decomposable_scoring",
        "uncalibrated_score_band",
        "phylogenetic_profile_cooccurrence",
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
            "uncalibrated_score_band": [False, False, False, False, True],
            "phylogenetic_profile_cooccurrence": [False, False, False, False, True],
            "matched_control_enrichment": [False, False, False, False, True],
            "archetype_clustering": [False, True, True, False, True],
            "structural_prioritization": [False, False, False, False, False],
            "fdr_controlled_candidates": [False, False, False, False, False],
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
        "T0_benchmark_validation_summary": summarize_benchmark_recovery(
            benchmark_recovery,
        ),
        "T1_benchmark_recovery": benchmark_recovery,
        "T2_top30_co_candidates": candidates.filter(pl.col("analyte") == "CO")
        .sort(
            "regulation_logit_score"
            if "regulation_logit_score" in candidates.columns
            else "candidate_score",
            descending=True,
        )
        .head(30),
        "T3_top30_cyd_control_candidates": candidates.filter(
            pl.col("analyte") == "cyd_control",
        )
        .sort(
            "regulation_logit_score"
            if "regulation_logit_score" in candidates.columns
            else "candidate_score",
            descending=True,
        )
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
