from __future__ import annotations

from pathlib import Path

import polars as pl

from gasregnet.reports.tables import tool_feature_comparison, write_publication_tables
from tests.unit.archetypes.test_cluster import _loci
from tests.unit.test_schemas import archetypes_frame, candidates_frame, enrichment_frame


def test_tool_feature_comparison_marks_gasregnet_features() -> None:
    matrix = tool_feature_comparison()
    gasregnet = matrix.filter(pl.col("tool") == "GasRegNet")

    assert gasregnet["matched_control_enrichment"].item() is True
    assert gasregnet["fdr_controlled_candidates"].item() is True
    assert gasregnet["operon_level_score_band"].item() is True
    assert gasregnet["phylogenetic_profile_cooccurrence"].item() is True


def test_write_publication_tables_outputs_csv_and_markdown(tmp_path: Path) -> None:
    benchmark = (
        _loci()
        .select(
            pl.col("locus_id").alias("benchmark_id"),
            "analyte",
        )
        .with_columns(
            pl.lit("CooA").alias("protein_name"),
            pl.lit("Org").alias("organism"),
            pl.lit(True).alias("hit"),
            pl.lit(1).alias("rank"),
            pl.lit(9.0).alias("candidate_score"),
        )
    )
    candidates = pl.concat(
        [
            candidates_frame(),
            candidates_frame().with_columns(
                pl.lit("CN").alias("analyte"),
                pl.lit("cand2").alias("candidate_id"),
            ),
        ],
    )

    outputs = write_publication_tables(
        benchmark_recovery=benchmark,
        candidates=candidates,
        enrichment=enrichment_frame(),
        archetypes=archetypes_frame(),
        out_dir=tmp_path,
    )

    assert set(outputs) == {
        "T1_benchmark_recovery",
        "T2_top30_co_candidates",
        "T3_top30_cyd_control_candidates",
        "T4_regulator_family_enrichment",
        "T5_archetype_catalog",
        "T6_tool_feature_comparison",
    }
    for csv_path, markdown_path in outputs.values():
        assert csv_path.exists()
        assert markdown_path.exists()
        assert markdown_path.read_text(encoding="utf-8").startswith("|")
    assert ",PAS," in outputs["T2_top30_co_candidates"][0].read_text(
        encoding="utf-8",
    )
