from __future__ import annotations

from pathlib import Path

import polars as pl

from gasregnet.reports.figures import (
    figure_1_workflow_and_recovery,
    figure_2_locus_landscape,
    figure_3_archetype_atlas,
    figure_4_chemistry_partition,
    figure_5_candidate_ranking,
    figure_6_structure_hypotheses,
)
from tests.unit.archetypes.test_cluster import _loci
from tests.unit.test_schemas import (
    archetypes_frame,
    candidates_frame,
    enrichment_frame,
)


def _assert_png_and_svg(png_path: Path) -> None:
    assert png_path.exists()
    assert png_path.stat().st_size > 0
    svg_path = png_path.with_suffix(".svg")
    assert svg_path.exists()
    assert svg_path.stat().st_size > 0


def test_publication_figure_builders_write_png_and_svg(tmp_path: Path) -> None:
    benchmark = pl.DataFrame(
        {
            "benchmark_id": ["co_known", "cn_known"],
            "analyte": ["CO", "CN"],
            "hit": [True, False],
        },
    )
    candidates = candidates_frame().with_columns(
        pl.lit(0.75).alias("structural_plausibility_score"),
    )

    paths = [
        figure_1_workflow_and_recovery(benchmark, tmp_path),
        figure_2_locus_landscape(_loci(), tmp_path),
        figure_3_archetype_atlas(archetypes_frame(), tmp_path),
        figure_4_chemistry_partition(enrichment_frame(), tmp_path),
        figure_5_candidate_ranking(candidates, tmp_path),
        figure_6_structure_hypotheses(candidates, tmp_path / "structures", tmp_path),
    ]

    for path in paths:
        _assert_png_and_svg(path)


def test_figure_builder_handles_empty_input(tmp_path: Path) -> None:
    path = figure_5_candidate_ranking(pl.DataFrame(), tmp_path)

    _assert_png_and_svg(path)
