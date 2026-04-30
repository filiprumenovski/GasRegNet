from __future__ import annotations

from pathlib import Path

import polars as pl

from scripts.run_sqlite_demo import run_sqlite_demo


def test_sqlite_demo_pipeline_writes_report_artifacts(tmp_path: Path) -> None:
    marker = run_sqlite_demo(
        out_dir=tmp_path / "repro",
        config_path=Path("configs/headline.yaml"),
        sqlite_path=tmp_path / "mini.sqlite",
    )

    out_dir = tmp_path / "repro"
    assert marker == out_dir / "README.txt"
    assert marker.exists()
    assert (out_dir / "manifest.json").exists()
    assert (out_dir / "config.resolved.yaml").exists()
    assert pl.read_parquet(out_dir / "intermediate" / "candidates.parquet").height == 2
    assert (out_dir / "tables" / "T6_tool_feature_comparison.md").exists()
    assert (out_dir / "figures" / "figure_4_chemistry_partition.png").exists()
    assert (
        out_dir / "captions" / "figure_1_workflow_and_recovery.md"
    ).read_text(encoding="utf-8").startswith("GasRegNet recovers")
