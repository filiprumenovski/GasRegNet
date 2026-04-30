from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest
from click.testing import CliRunner

from gasregnet.cli import app
from scripts.run_sqlite_demo import run_sqlite_demo

runner = CliRunner()


def test_help_lists_core_subcommands() -> None:
    result = runner.invoke(app, ["--help"])

    assert result.exit_code == 0
    for command in (
        "validate-config",
        "build-benchmark",
        "run-sqlite",
        "diamond-search",
        "annotate",
        "score",
        "enrich",
        "archetypes",
        "report",
        "repro",
    ):
        assert command in result.output


def test_validate_config_writes_run_metadata(tmp_path: Path) -> None:
    result = runner.invoke(
        app,
        [
            "validate-config",
            "--config",
            "configs",
            "--out",
            str(tmp_path),
        ],
    )

    assert result.exit_code == 0
    assert "config valid" in result.output
    assert (tmp_path / "config.resolved.yaml").exists()
    assert (tmp_path / "manifest.json").exists()


def test_build_benchmark_writes_csv_header(tmp_path: Path) -> None:
    out = tmp_path / "benchmark.csv"

    result = runner.invoke(app, ["build-benchmark", "--out", str(out)])

    assert result.exit_code == 0
    assert out.read_text(encoding="utf-8").startswith("benchmark_id,analyte")


def test_report_command_writes_artifacts(tmp_path: Path) -> None:
    results = tmp_path / "results"
    report = tmp_path / "report"
    run_sqlite_demo(
        out_dir=results,
        config_path=Path("configs/headline.yaml"),
        sqlite_path=tmp_path / "mini.sqlite",
    )

    result = runner.invoke(
        app,
        [
            "report",
            "--results",
            str(results),
            "--out",
            str(report),
        ],
    )

    assert result.exit_code == 0
    assert (report / "tables" / "T6_tool_feature_comparison.csv").exists()
    assert (report / "figures" / "figure_5_candidate_ranking.png").exists()
    assert (report / "captions" / "figure_4_chemistry_partition.md").exists()


def test_repro_command_invokes_snakemake(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls: list[list[str]] = []

    def fake_run(command: list[str], **_: Any) -> None:
        calls.append(command)

    monkeypatch.setattr("subprocess.run", fake_run)

    result = runner.invoke(
        app,
        [
            "repro",
            "--config",
            "configs/headline.yaml",
            "--out",
            str(tmp_path / "headline"),
            "--sqlite",
            "tests/fixtures/mini_efi.sqlite",
        ],
    )

    assert result.exit_code == 0
    assert calls
    assert calls[0][:5] == [
        "uv",
        "run",
        "snakemake",
        "-s",
        "workflows/sqlite_mode.smk",
    ]
    assert f"out_dir={tmp_path / 'headline'}" in calls[0]


def test_placeholder_command_exits_clearly() -> None:
    result = runner.invoke(app, ["annotate"])

    assert result.exit_code == 2
    assert "not implemented" in result.output
