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
        "fetch-assets",
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


def test_fetch_assets_command_writes_manifest_outputs(tmp_path: Path) -> None:
    source = tmp_path / "source.faa"
    source.write_text(">seed\nMA\n", encoding="utf-8")
    manifest = tmp_path / "assets.yaml"
    manifest.write_text(
        f"""
assets:
  - name: seed
    output: imported/seed.faa
    urls:
      - {source.as_uri()}
""",
        encoding="utf-8",
    )

    result = runner.invoke(
        app,
        [
            "fetch-assets",
            "--manifest",
            str(manifest),
            "--root",
            str(tmp_path),
            "--downloader",
            "urllib",
            "--force",
        ],
    )

    assert result.exit_code == 0
    assert (tmp_path / "imported" / "seed.faa").read_text(encoding="utf-8") == (
        ">seed\nMA\n"
    )


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


def test_stage_commands_write_scored_enrichment_and_archetypes(tmp_path: Path) -> None:
    source = tmp_path / "source"
    annotated = tmp_path / "annotated"
    scored = tmp_path / "scored"
    enriched = tmp_path / "enriched"
    archetypes = tmp_path / "archetypes"
    run_sqlite_demo(
        out_dir=source,
        config_path=Path("configs/headline.yaml"),
        sqlite_path=tmp_path / "mini.sqlite",
    )

    annotate_result = runner.invoke(
        app,
        [
            "annotate",
            "--neighborhoods",
            str(source),
            "--config",
            "configs",
            "--out",
            str(annotated),
        ],
    )
    score_result = runner.invoke(
        app,
        [
            "score",
            "--neighborhoods",
            str(annotated),
            "--config",
            "configs",
            "--out",
            str(scored),
        ],
    )
    enrich_result = runner.invoke(
        app,
        [
            "enrich",
            "--scored",
            str(scored),
            "--config",
            "configs",
            "--out",
            str(enriched),
        ],
    )
    archetypes_result = runner.invoke(
        app,
        [
            "archetypes",
            "--scored",
            str(scored),
            "--config",
            "configs",
            "--out",
            str(archetypes),
        ],
    )

    assert annotate_result.exit_code == 0
    assert score_result.exit_code == 0
    assert enrich_result.exit_code == 0
    assert archetypes_result.exit_code == 0
    assert (annotated / "intermediate" / "genes.parquet").exists()
    assert (scored / "intermediate" / "candidates.parquet").exists()
    assert (enriched / "intermediate" / "enrichment.parquet").exists()
    assert (archetypes / "intermediate" / "archetypes.parquet").exists()


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


def test_annotate_command_accepts_optional_annotation_tables(tmp_path: Path) -> None:
    source = tmp_path / "source"
    out = tmp_path / "annotated"
    pfam = tmp_path / "pfam.csv"
    interpro = tmp_path / "interpro.csv"
    run_sqlite_demo(
        out_dir=source,
        config_path=Path("configs/headline.yaml"),
        sqlite_path=tmp_path / "mini.sqlite",
    )
    pfam.write_text(
        "gene_accession,pfam_id,pfam_description\nREG_UP,PF01047,HTH\n",
        encoding="utf-8",
    )
    interpro.write_text(
        "gene_accession,interpro_id,interpro_description\n",
        encoding="utf-8",
    )

    result = runner.invoke(
        app,
        [
            "annotate",
            "--neighborhoods",
            str(source),
            "--config",
            "configs",
            "--out",
            str(out),
            "--pfam-table",
            str(pfam),
            "--interpro-table",
            str(interpro),
        ],
    )

    assert result.exit_code == 0
    assert (out / "intermediate" / "genes.parquet").exists()
