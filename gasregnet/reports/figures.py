"""Publication figure builders."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import polars as pl

STYLE_PATH = Path(__file__).with_name("gasregnet.mplstyle")


def _save(fig: Any, out_dir: Path, stem: str) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    png_path = out_dir / f"{stem}.png"
    svg_path = out_dir / f"{stem}.svg"
    fig.tight_layout()
    fig.savefig(png_path, dpi=300)
    fig.savefig(svg_path)
    plt.close(fig)
    return png_path


def _empty_figure(out_dir: Path, stem: str, title: str) -> Path:
    with plt.style.context(STYLE_PATH):
        fig, ax = plt.subplots(figsize=(5, 3))
        ax.set_axis_off()
        ax.text(0.5, 0.5, "No data", ha="center", va="center")
        ax.set_title(title, loc="left")
        return _save(fig, out_dir, stem)


def figure_1_workflow_and_recovery(
    benchmark_results: pl.DataFrame,
    out_dir: Path,
) -> Path:
    """Render benchmark recovery by analyte."""

    if benchmark_results.is_empty() or "hit" not in benchmark_results.columns:
        return _empty_figure(
            out_dir,
            "figure_1_workflow_and_recovery",
            "Benchmark recovery",
        )
    summary = (
        benchmark_results.group_by("analyte")
        .agg(pl.col("hit").mean().mul(100).alias("recall_percent"))
        .sort("analyte")
    )
    with plt.style.context(STYLE_PATH):
        fig, ax = plt.subplots(figsize=(5.5, 3.2))
        ax.bar(summary["analyte"].to_list(), summary["recall_percent"].to_list())
        ax.set_ylim(0, 100)
        ax.set_ylabel("Recall (%)")
        ax.set_title("Benchmark sensor recovery", loc="left")
        return _save(fig, out_dir, "figure_1_workflow_and_recovery")


def figure_2_locus_landscape(loci: pl.DataFrame, out_dir: Path) -> Path:
    """Render locus-score distributions by analyte and confidence."""

    if loci.is_empty():
        return _empty_figure(out_dir, "figure_2_locus_landscape", "Locus landscape")
    colors = {"high": "#4C78A8", "medium": "#F58518", "low": "#BAB0AC"}
    with plt.style.context(STYLE_PATH):
        fig, ax = plt.subplots(figsize=(6.2, 3.6))
        sorted_loci = loci.sort(["analyte", "locus_score"])
        for index, row in enumerate(sorted_loci.iter_rows(named=True)):
            confidence = str(row.get("locus_confidence", "low"))
            ax.scatter(
                index,
                float(row["locus_score"]),
                color=colors.get(confidence, "#BAB0AC"),
                s=42,
            )
        ax.set_ylabel("Locus score")
        ax.set_xlabel("Loci sorted by analyte and score")
        ax.set_title("CO and HCN locus landscape", loc="left")
        return _save(fig, out_dir, "figure_2_locus_landscape")


def figure_3_archetype_atlas(archetypes: pl.DataFrame, out_dir: Path) -> Path:
    """Render recurrent archetype abundance."""

    if archetypes.is_empty():
        return _empty_figure(out_dir, "figure_3_archetype_atlas", "Archetype atlas")
    top = archetypes.sort("n_loci", descending=True).head(12)
    labels = [
        f"{row['archetype_id']} ({row['analyte']})"
        for row in top.iter_rows(named=True)
    ]
    with plt.style.context(STYLE_PATH):
        fig, ax = plt.subplots(figsize=(7, 4))
        ax.barh(labels[::-1], top["n_loci"].to_list()[::-1])
        ax.set_xlabel("Loci")
        ax.set_title("Recurrent neighborhood archetypes", loc="left")
        return _save(fig, out_dir, "figure_3_archetype_atlas")


def figure_4_chemistry_partition(enrichment: pl.DataFrame, out_dir: Path) -> Path:
    """Render enrichment effect sizes with q-value scaling."""

    if enrichment.is_empty():
        return _empty_figure(
            out_dir,
            "figure_4_chemistry_partition",
            "Chemistry partition",
        )
    with plt.style.context(STYLE_PATH):
        fig, ax = plt.subplots(figsize=(6.5, 3.8))
        for row in enrichment.sort("q_value").iter_rows(named=True):
            label = f"{row['analyte']} {row['feature_name']}"
            ax.scatter(
                float(row["odds_ratio"]),
                label,
                s=max(30.0, 180.0 * (1.0 - float(row["q_value"]))),
            )
        ax.set_xlabel("Odds ratio")
        ax.set_title("Regulator feature enrichment", loc="left")
        return _save(fig, out_dir, "figure_4_chemistry_partition")


def figure_5_candidate_ranking(candidates: pl.DataFrame, out_dir: Path) -> Path:
    """Render top candidate sensor scores."""

    if candidates.is_empty():
        return _empty_figure(out_dir, "figure_5_candidate_ranking", "Candidate ranking")
    top = candidates.sort("candidate_score", descending=True).head(30)
    with plt.style.context(STYLE_PATH):
        fig, ax = plt.subplots(figsize=(7, 4.5))
        ax.barh(
            top["candidate_id"].to_list()[::-1],
            top["candidate_score"].to_list()[::-1],
        )
        ax.set_xlabel("Candidate score")
        ax.set_title("Top candidate gas sensors", loc="left")
        return _save(fig, out_dir, "figure_5_candidate_ranking")


def figure_6_structure_hypotheses(
    top_candidates: pl.DataFrame,
    structures_dir: Path,
    out_dir: Path,
) -> Path:
    """Render structure-prioritized candidates."""

    del structures_dir
    if (
        top_candidates.is_empty()
        or "structural_plausibility_score" not in top_candidates.columns
    ):
        return _empty_figure(
            out_dir,
            "figure_6_structure_hypotheses",
            "Structure hypotheses",
        )
    scored = top_candidates.filter(
        pl.col("structural_plausibility_score").is_not_null(),
    )
    if scored.is_empty():
        return _empty_figure(
            out_dir,
            "figure_6_structure_hypotheses",
            "Structure hypotheses",
        )
    top = scored.sort("structural_plausibility_score", descending=True).head(10)
    with plt.style.context(STYLE_PATH):
        fig, ax = plt.subplots(figsize=(6.5, 3.8))
        ax.barh(
            top["candidate_id"].to_list()[::-1],
            top["structural_plausibility_score"].to_list()[::-1],
        )
        ax.set_xlabel("Structural plausibility")
        ax.set_title("Structure-guided sensor hypotheses", loc="left")
        return _save(fig, out_dir, "figure_6_structure_hypotheses")
