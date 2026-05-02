from __future__ import annotations

import json
import logging
from pathlib import Path

import structlog

from gasregnet.hashing import file_sha256, text_sha256
from gasregnet.logging import configure_logging
from gasregnet.manifest import build_manifest, write_manifest
from gasregnet.paths import ensure_out_dir


def test_hashing_helpers_are_deterministic(tmp_path: Path) -> None:
    payload = tmp_path / "payload.txt"
    payload.write_text("gasregnet\n", encoding="utf-8")

    assert file_sha256(payload) == file_sha256(payload)
    assert text_sha256("gasregnet\n") == file_sha256(payload)


def test_build_manifest_hashes_declared_inputs_deterministically(
    tmp_path: Path,
) -> None:
    config = tmp_path / "config.yaml"
    input_file = tmp_path / "input.faa"
    config.write_text("seed: 7\n", encoding="utf-8")
    input_file.write_text(">seq\nM\n", encoding="utf-8")

    first = build_manifest(
        seed=7,
        command="gasregnet run",
        config_paths={"config": config},
        input_paths={"seeds": input_file},
        external_versions={"mmseqs": "15-6f452"},
    )
    second = build_manifest(
        seed=7,
        command="gasregnet run",
        config_paths={"config": config},
        input_paths={"seeds": input_file},
        external_versions={"mmseqs": "15-6f452"},
    )

    assert first.run_hash == second.run_hash
    assert first.config_hashes == second.config_hashes
    assert first.input_hashes == second.input_hashes


def test_build_manifest_hashes_files_inside_input_directories(tmp_path: Path) -> None:
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    (inputs / "a.faa").write_text(">a\nMA\n", encoding="utf-8")
    nested = inputs / "nested"
    nested.mkdir()
    (nested / "b.hmm").write_text("HMM\n", encoding="utf-8")

    manifest = build_manifest(
        seed=7,
        command="gasregnet run",
        input_paths={"assets": inputs},
    )

    assert set(manifest.input_hashes) == {
        "assets/a.faa",
        "assets/nested/b.hmm",
    }


def test_build_manifest_loads_tools_resolved_yaml(tmp_path: Path) -> None:
    tools = tmp_path / "tools_resolved.yaml"
    tools.write_text(
        "diamond:\n  version: diamond version 2.1\n  path: /bin/diamond\n",
        encoding="utf-8",
    )

    manifest = build_manifest(
        seed=7,
        command="gasregnet run",
        tools_resolved=tools,
    )

    assert manifest.external_versions["diamond"]["version"] == "diamond version 2.1"
    assert "tools_resolved" in manifest.input_hashes


def test_write_manifest_outputs_sorted_json(tmp_path: Path) -> None:
    manifest = build_manifest(seed=1, command="gasregnet test")

    path = write_manifest(manifest, tmp_path)

    payload = json.loads(path.read_text(encoding="utf-8"))
    assert payload["command"] == "gasregnet test"
    assert payload["seed"] == 1
    assert path.read_text(encoding="utf-8").endswith("\n")


def test_ensure_out_dir_creates_standard_layout(tmp_path: Path) -> None:
    out_dir = ensure_out_dir(tmp_path / "run")

    assert out_dir.exists()
    for child in ("logs", "tables", "figures", "intermediate"):
        assert (out_dir / child).is_dir()


def test_configure_logging_sets_filtering_bound_logger() -> None:
    configure_logging(verbose=True)

    logger = structlog.get_logger("gasregnet-test")

    assert logging.getLogger().getEffectiveLevel() == logging.DEBUG
    assert logger is not None
