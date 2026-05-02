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
        "check-tools",
        "fetch-assets",
        "build-profiles",
        "build-seed-databases",
        "index-refseq",
        "index-refseq-corpus",
        "enumerate-shards",
        "query-refseq",
        "query-refseq-corpus",
        "summarize-refseq-corpus",
        "scan-refseq-corpus",
        "detect-anchors",
        "detect-anchors-profile",
        "extract-neighborhoods",
        "evaluate-benchmark",
        "build-benchmark",
        "run-sqlite",
        "diamond-search",
        "annotate",
        "assign-roles",
        "score",
        "simulate-synthetic-truth",
        "enrich",
        "archetypes",
        "report",
        "repro",
        "corpus-discovery",
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


def test_check_tools_command_reports_missing_required_tools(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    monkeypatch.setenv("PATH", "")

    result = runner.invoke(app, ["check-tools", "--out", str(tmp_path / "tools.yaml")])

    assert result.exit_code != 0
    assert "required external tool not found" in result.output


def test_build_seed_databases_runs_diamond_per_family(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    calls: list[list[str]] = []

    def fake_which(name: str) -> str | None:
        return "/usr/local/bin/diamond" if name == "diamond" else None

    def fake_run(
        args: list[str],
        *,
        check: bool,
        capture_output: bool,
        text: bool,
    ) -> None:
        calls.append(args)
        assert check is True
        assert capture_output is True
        assert text is True

    monkeypatch.setattr("gasregnet.cli.shutil.which", fake_which)
    monkeypatch.setattr("gasregnet.cli.subprocess.run", fake_run)

    result = runner.invoke(
        app,
        [
            "build-seed-databases",
            "--config",
            "configs/headline.yaml",
            "--out",
            str(tmp_path),
        ],
    )

    assert result.exit_code == 0
    assert len(calls) == 19
    assert [
        "/usr/local/bin/diamond",
        "makedb",
        "--in",
        "data/seeds/co_coxL_anchor_seeds.faa",
        "--db",
        str(tmp_path / "CO__coxL.dmnd"),
    ] in calls
    assert str(tmp_path / "cyd_control__cydA.dmnd") in result.output


def test_build_benchmark_writes_csv_header(tmp_path: Path) -> None:
    out = tmp_path / "benchmark.csv"

    result = runner.invoke(app, ["build-benchmark", "--out", str(out)])

    assert result.exit_code == 0
    assert out.read_text(encoding="utf-8").startswith("benchmark_id,analyte")


def test_build_benchmark_v2_writes_regulator_benchmark(tmp_path: Path) -> None:
    out = tmp_path / "regulators_v2.csv"

    result = runner.invoke(
        app,
        ["build-benchmark", "--version", "v2", "--out", str(out)],
    )

    assert result.exit_code == 0
    text = out.read_text(encoding="utf-8")
    assert "CooA" in text
    assert "verify_pmid" in text


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
    roles = tmp_path / "roles"
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
    roles_result = runner.invoke(
        app,
        [
            "assign-roles",
            "--neighborhoods",
            str(annotated),
            "--config",
            "configs",
            "--out",
            str(roles),
        ],
    )

    assert annotate_result.exit_code == 0
    assert score_result.exit_code == 0
    assert enrich_result.exit_code == 0
    assert archetypes_result.exit_code == 0
    assert roles_result.exit_code == 0
    assert (annotated / "intermediate" / "genes.parquet").exists()
    assert (scored / "intermediate" / "candidates.parquet").exists()
    assert (enriched / "intermediate" / "enrichment.parquet").exists()
    assert (archetypes / "intermediate" / "archetypes.parquet").exists()
    assert (roles / "intermediate" / "sensor_regulator_pairs.parquet").exists()


def test_simulate_synthetic_truth_command_writes_truth_frames(tmp_path: Path) -> None:
    result = runner.invoke(
        app,
        [
            "simulate-synthetic-truth",
            "--out",
            str(tmp_path),
            "--n-genomes",
            "8",
            "--annotation-noise",
            "0",
        ],
    )

    assert result.exit_code == 0
    assert (tmp_path / "intermediate" / "loci.parquet").exists()
    assert (tmp_path / "intermediate" / "genes.parquet").exists()
    assert (tmp_path / "intermediate" / "synthetic_ground_truth.csv").exists()


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


def test_corpus_discovery_command_invokes_snakemake(
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
            "corpus-discovery",
            "--corpus-config",
            "configs/corpus_discovery.yaml",
            "--config",
            "configs/headline.yaml",
            "--out",
            str(tmp_path / "corpus"),
        ],
    )

    assert result.exit_code == 0
    assert calls
    assert calls[0][:5] == [
        "uv",
        "run",
        "snakemake",
        "-s",
        "workflows/corpus_discovery.smk",
    ]
    assert f"out_dir={tmp_path / 'corpus'}" in calls[0]


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
