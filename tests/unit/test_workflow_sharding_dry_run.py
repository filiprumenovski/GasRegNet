from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import polars as pl


def _write_store(store: Path) -> None:
    store.mkdir(parents=True)
    pl.DataFrame(
        {
            "dataset_name": [
                "actino_1",
                "bacillo_1",
                "pseudo_1",
                "pseudo_2",
            ],
            "phylum": [
                "Actinomycetota",
                "Bacillota",
                "Pseudomonadota",
                "Pseudomonadota",
            ],
            "n_proteins": [10, 20, 30, 40],
        },
    ).write_parquet(store / "datasets.parquet")


def test_corpus_workflow_dry_run_fans_out_by_phylum(tmp_path: Path) -> None:
    repo = Path(__file__).parents[2]
    store = tmp_path / "store"
    _write_store(store)

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "snakemake",
            "--snakefile",
            "workflows/corpus_discovery.smk",
            "--dry-run",
            "--profile",
            "workflows/profiles/local",
            "--config",
            f"out_dir={tmp_path / 'out'}",
            "corpus_config=configs/corpus_discovery.yaml",
            "config_path=configs/headline.yaml",
            "root=.",
            f"store={store}",
            "sharding_strategy=by_phylum",
        ],
        cwd=repo,
        capture_output=True,
        check=False,
        text=True,
    )

    assert result.returncode == 0, result.stderr
    assert (store / "shards.parquet").exists()
    assert result.stdout.count("rule detect_anchors_shard:") == 3


def test_corpus_workflow_dry_run_fans_out_by_n_genomes(tmp_path: Path) -> None:
    repo = Path(__file__).parents[2]
    store = tmp_path / "store"
    _write_store(store)

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "snakemake",
            "--snakefile",
            "workflows/corpus_discovery.smk",
            "--dry-run",
            "--cores",
            "1",
            "--config",
            f"out_dir={tmp_path / 'out'}",
            "corpus_config=configs/corpus_discovery.yaml",
            "config_path=configs/headline.yaml",
            "root=.",
            f"store={store}",
            "sharding_strategy=by_n_genomes",
            "n_genomes_per_shard=2",
        ],
        cwd=repo,
        capture_output=True,
        check=False,
        text=True,
    )

    assert result.returncode == 0, result.stderr
    assert result.stdout.count("rule detect_anchors_shard:") == 2
