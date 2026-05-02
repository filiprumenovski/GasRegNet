from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

WORKFLOWS = [
    pytest.param(
        "workflows/sqlite_mode.smk",
        {
            "out_dir": "sqlite",
            "config_path": "configs/headline.yaml",
            "sqlite": "tests/fixtures/mini_efi.sqlite",
        },
        id="sqlite",
    ),
    pytest.param(
        "workflows/full_discovery.smk",
        {
            "out_dir": "full",
            "config_path": "configs/headline.yaml",
            "sqlite": "tests/fixtures/mini_efi.sqlite",
        },
        id="full-discovery",
    ),
    pytest.param(
        "workflows/corpus_discovery.smk",
        {
            "out_dir": "corpus",
            "corpus_config": "configs/corpus_discovery.yaml",
            "config_path": "configs/headline.yaml",
            "root": ".",
        },
        id="corpus-discovery",
    ),
    pytest.param(
        "workflows/diamond_mode.smk",
        {
            "out_dir": "diamond",
            "query": "data/seeds/co_anchor_seeds.faa",
            "db": "tests/fixtures/mini_efi.sqlite",
        },
        id="diamond",
    ),
]


@pytest.mark.parametrize(("workflow", "config"), WORKFLOWS)
def test_snakemake_workflow_dry_run(
    tmp_path: Path,
    workflow: str,
    config: dict[str, str],
) -> None:
    resolved_config = {
        key: str(tmp_path / value) if key == "out_dir" else value
        for key, value in config.items()
    }
    config_args = [f"{key}={value}" for key, value in resolved_config.items()]

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "snakemake",
            "--snakefile",
            workflow,
            "--dry-run",
            "--cores",
            "1",
            "--config",
            *config_args,
        ],
        cwd=Path(__file__).parents[2],
        capture_output=True,
        check=False,
        text=True,
    )

    assert result.returncode == 0, result.stderr
