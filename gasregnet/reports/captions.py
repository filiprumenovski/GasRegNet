"""Result-led figure captions."""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path
from typing import Any

import polars as pl


def _unique_count(frame: pl.DataFrame, column: str) -> int:
    if column not in frame.columns or frame.is_empty():
        return 0
    return int(frame.select(pl.col(column).n_unique()).item())


def _count_value(frame: pl.DataFrame, column: str, value: object) -> int:
    if column not in frame.columns or frame.is_empty():
        return 0
    return frame.filter(pl.col(column) == value).height


def _percent(numerator: float, denominator: float) -> str:
    if denominator == 0:
        return "0.0%"
    return f"{(numerator / denominator) * 100:.1f}%"


def _top_row(frame: pl.DataFrame, sort_by: list[str]) -> dict[str, Any] | None:
    if frame.is_empty():
        return None
    return frame.sort(sort_by).row(0, named=True)


def _float_text(value: Any, digits: int = 2) -> str:
    if value is None:
        return "not scored"
    return f"{float(value):.{digits}g}"


def figure_1_workflow_and_recovery_caption(benchmark_results: pl.DataFrame) -> str:
    """Caption for benchmark recovery and workflow figure."""

    recovered = _count_value(benchmark_results, "hit", True)
    total = benchmark_results.height
    organisms = _unique_count(benchmark_results, "organism")
    return (
        "GasRegNet recovers "
        f"{recovered} of {total} benchmark bacterial gas sensors across "
        f"{organisms} organisms (recall {_percent(recovered, total)})."
    )


def figure_2_locus_landscape_caption(loci: pl.DataFrame) -> str:
    """Caption for CO/CN locus landscape figure."""

    co_loci = _count_value(loci, "analyte", "CO")
    cn_loci = _count_value(loci, "analyte", "CN")
    high_confidence = _count_value(loci, "locus_confidence", "high")
    taxa = _unique_count(loci, "taxon_id")
    return (
        "CO and HCN searches resolve "
        f"{co_loci} CO loci and {cn_loci} HCN loci across {taxa} taxa, "
        f"with {high_confidence} high-confidence neighborhoods."
    )


def figure_3_archetype_atlas_caption(archetypes: pl.DataFrame) -> str:
    """Caption for recurrent gene-architecture atlas."""

    clauses: list[str] = []
    for analyte in ("CO", "CN"):
        subset = archetypes.filter(pl.col("analyte") == analyte)
        total_loci = float(subset.select(pl.col("n_loci").sum()).item() or 0)
        top_loci = float(
            subset.sort("n_loci", descending=True)
            .head(6)
            .select(pl.col("n_loci").sum())
            .item()
            or 0,
        )
        clauses.append(
            f"{analyte}: {subset.height} archetypes covering "
            f"{_percent(top_loci, total_loci)} of assigned loci",
        )
    return (
        "Recurrent gene-neighborhood archetypes are compact: "
        + "; ".join(clauses)
        + "."
    )


def figure_4_chemistry_partition_caption(enrichment: pl.DataFrame) -> str:
    """Caption for regulator-family and sensory-domain enrichment figure."""

    top = _top_row(enrichment, ["q_value", "p_value"])
    if top is None:
        return "Sensory chemistry partitioning has no enrichment rows in this run."
    return (
        "Sensory chemistry partitions regulator features, led by "
        f"{top['analyte']} {top['feature_type']} {top['feature_name']} "
        f"(odds ratio {_float_text(top['odds_ratio'])}, "
        f"q={_float_text(top['q_value'])})."
    )


def figure_5_candidate_ranking_caption(candidates: pl.DataFrame) -> str:
    """Caption for decomposable candidate ranking figure."""

    if "candidate_score_q" in candidates.columns:
        high_confidence = candidates.filter(pl.col("candidate_score_q") <= 0.05).height
    else:
        high_confidence = 0
    if candidates.is_empty():
        return "GasRegNet nominated 0 high-confidence candidate sensors in this run."
    top = candidates.sort("candidate_score", descending=True).row(0, named=True)
    return (
        "GasRegNet nominates "
        f"{high_confidence} high-confidence candidate sensors, headed by "
        f"{top['candidate_id']} in {top['organism']} "
        f"(score {_float_text(top['candidate_score'])})."
    )


def figure_6_structure_hypotheses_caption(top_candidates: pl.DataFrame) -> str:
    """Caption for structure-guided sensor hypotheses figure."""

    if "structural_plausibility_score" not in top_candidates.columns:
        return (
            "Structural prioritization has 0 candidates with structure-derived "
            "scores."
        )
    scored = top_candidates.filter(
        pl.col("structural_plausibility_score").is_not_null(),
    )
    if scored.is_empty():
        return (
            "Structural prioritization has 0 candidates with structure-derived "
            "scores."
        )
    top_dict = scored.sort("structural_plausibility_score", descending=True).row(
        0,
        named=True,
    )
    return (
        "Structural prioritization supports "
        f"{scored.height} candidate sensor hypotheses, led by "
        f"{top_dict['candidate_id']} "
        f"(structure score {_float_text(top_dict['structural_plausibility_score'])})."
    )


def build_result_led_captions(
    *,
    benchmark_results: pl.DataFrame,
    loci: pl.DataFrame,
    archetypes: pl.DataFrame,
    enrichment: pl.DataFrame,
    candidates: pl.DataFrame,
    top_candidates: pl.DataFrame,
) -> dict[str, str]:
    """Build all main-figure captions from run outputs."""

    return {
        "figure_1_workflow_and_recovery": figure_1_workflow_and_recovery_caption(
            benchmark_results,
        ),
        "figure_2_locus_landscape": figure_2_locus_landscape_caption(loci),
        "figure_3_archetype_atlas": figure_3_archetype_atlas_caption(archetypes),
        "figure_4_chemistry_partition": figure_4_chemistry_partition_caption(
            enrichment,
        ),
        "figure_5_candidate_ranking": figure_5_candidate_ranking_caption(candidates),
        "figure_6_structure_hypotheses": figure_6_structure_hypotheses_caption(
            top_candidates,
        ),
    }


def write_caption_files(captions: Mapping[str, str], out_dir: Path) -> dict[str, Path]:
    """Write one Markdown caption file per figure."""

    out_dir.mkdir(parents=True, exist_ok=True)
    outputs: dict[str, Path] = {}
    for stem, caption in captions.items():
        path = out_dir / f"{stem}.md"
        path.write_text(caption + "\n", encoding="utf-8")
        outputs[stem] = path
    return outputs
