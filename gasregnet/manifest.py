"""Run manifest creation."""

from __future__ import annotations

import json
import platform
from dataclasses import asdict, dataclass, field
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from gasregnet import __version__
from gasregnet.check_tools import read_tools_resolved
from gasregnet.hashing import file_sha256, text_sha256


@dataclass(frozen=True)
class RunManifest:
    run_hash: str
    seed: int
    command: str
    package_version: str = __version__
    timestamp: str = field(default_factory=lambda: datetime.now(UTC).isoformat())
    python_version: str = field(default_factory=platform.python_version)
    config_hashes: dict[str, str] = field(default_factory=dict)
    input_hashes: dict[str, str] = field(default_factory=dict)
    external_versions: dict[str, Any] = field(default_factory=dict)
    wall_clock_seconds: float | None = None


def _hash_paths(paths: dict[str, Path]) -> dict[str, str]:
    hashes: dict[str, str] = {}
    for name, path in sorted(paths.items()):
        if path.is_dir():
            children = sorted(
                candidate for candidate in path.rglob("*") if candidate.is_file()
            )
            for child in children:
                relative = child.relative_to(path).as_posix()
                hashes[f"{name}/{relative}"] = file_sha256(child)
            continue
        hashes[name] = file_sha256(path)
    return hashes


def build_manifest(
    *,
    seed: int,
    command: str,
    config_paths: dict[str, Path] | None = None,
    input_paths: dict[str, Path] | None = None,
    external_versions: dict[str, Any] | None = None,
    tools_resolved: Path | None = None,
    wall_clock_seconds: float | None = None,
) -> RunManifest:
    """Build a run manifest with deterministic input/config hashes."""

    versions: dict[str, Any] = dict(sorted((external_versions or {}).items()))
    if tools_resolved is not None:
        versions.update(read_tools_resolved(tools_resolved))
        input_paths = {**(input_paths or {}), "tools_resolved": tools_resolved}
    config_hashes = _hash_paths(config_paths or {})
    input_hashes = _hash_paths(input_paths or {})
    run_hash = text_sha256(
        json.dumps(
            {
                "command": command,
                "config_hashes": config_hashes,
                "external_versions": versions,
                "input_hashes": input_hashes,
                "package_version": __version__,
                "seed": seed,
            },
            sort_keys=True,
        ),
    )
    return RunManifest(
        run_hash=run_hash,
        seed=seed,
        command=command,
        config_hashes=config_hashes,
        input_hashes=input_hashes,
        external_versions=versions,
        wall_clock_seconds=wall_clock_seconds,
    )


def write_manifest(manifest: RunManifest, out_dir: Path) -> Path:
    """Write a manifest to ``<out>/manifest.json``."""

    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "manifest.json"
    payload: dict[str, Any] = asdict(manifest)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return path
