from __future__ import annotations

from pathlib import Path

import polars as pl

from gasregnet.reports.captions import (
    build_result_led_captions,
    figure_1_workflow_and_recovery_caption,
    figure_5_candidate_ranking_caption,
    write_caption_files,
)
from tests.unit.archetypes.test_cluster import _loci
from tests.unit.test_schemas import (
    archetypes_frame,
    candidates_frame,
    enrichment_frame,
)


def test_figure_1_caption_uses_benchmark_counts() -> None:
    benchmark_results = pl.DataFrame(
        {
            "benchmark_id": ["known_a", "known_b"],
            "organism": ["Org1", "Org2"],
            "hit": [True, False],
        },
    )

    caption = figure_1_workflow_and_recovery_caption(benchmark_results)

    assert "1 of 2" in caption
    assert "2 organisms" in caption
    assert "50.0%" in caption


def test_figure_5_caption_reports_top_candidate() -> None:
    candidates = pl.concat(
        [
            candidates_frame().with_columns(pl.lit(10.0).alias("candidate_score")),
            candidates_frame().with_columns(
                pl.lit("cand2").alias("candidate_id"),
                pl.lit(20.0).alias("candidate_score"),
                pl.lit(0.20).alias("candidate_score_q"),
            ),
        ],
    )

    caption = figure_5_candidate_ranking_caption(candidates)

    assert "1 high-confidence" in caption
    assert "cand2" in caption
    assert "score 20" in caption


def test_build_and_write_all_caption_files(tmp_path: Path) -> None:
    candidates = candidates_frame().with_columns(
        pl.lit(0.8).alias("structural_plausibility_score"),
    )

    captions = build_result_led_captions(
        benchmark_results=pl.DataFrame(
            {
                "benchmark_id": ["known_a"],
                "organism": ["Org1"],
                "hit": [True],
            },
        ),
        loci=_loci(),
        archetypes=archetypes_frame(),
        enrichment=enrichment_frame(),
        candidates=candidates,
        top_candidates=candidates,
    )
    outputs = write_caption_files(captions, tmp_path)

    assert set(outputs) == {
        "figure_1_workflow_and_recovery",
        "figure_2_locus_landscape",
        "figure_3_archetype_atlas",
        "figure_4_chemistry_partition",
        "figure_5_candidate_ranking",
        "figure_6_structure_hypotheses",
    }
    assert (
        outputs["figure_4_chemistry_partition"]
        .read_text(
            encoding="utf-8",
        )
        .startswith("Sensory chemistry partitions")
    )
