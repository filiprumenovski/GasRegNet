"""Run manifest creation."""

from __future__ import annotations

import json
import platform
from dataclasses import asdict, dataclass, field
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from gasregnet import __version__


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
    external_versions: dict[str, str] = field(default_factory=dict)
    wall_clock_seconds: float | None = None


def write_manifest(manifest: RunManifest, out_dir: Path) -> Path:
    """Write a manifest to ``<out>/manifest.json``."""

    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "manifest.json"
    payload: dict[str, Any] = asdict(manifest)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return path
