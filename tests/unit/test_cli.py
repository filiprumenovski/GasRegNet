from __future__ import annotations

from pathlib import Path

from click.testing import CliRunner

from gasregnet.cli import app

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


def test_placeholder_command_exits_clearly() -> None:
    result = runner.invoke(app, ["report"])

    assert result.exit_code == 2
    assert "not implemented" in result.output
