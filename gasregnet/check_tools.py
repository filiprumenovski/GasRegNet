"""External binary provenance checks."""

from __future__ import annotations

import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml  # type: ignore[import-untyped]

from gasregnet.hashing import file_sha256


@dataclass(frozen=True)
class ToolProbe:
    name: str
    command: tuple[str, ...]
    required: bool = True


DEFAULT_TOOL_PROBES: tuple[ToolProbe, ...] = (
    ToolProbe("hmmsearch", ("hmmsearch", "-h")),
    ToolProbe("diamond", ("diamond", "version")),
    ToolProbe("foldseek", ("foldseek", "version"), required=False),
    ToolProbe("mafft", ("mafft", "--version"), required=False),
    ToolProbe("meme", ("meme", "-version"), required=False),
)


def _first_version_line(text: str) -> str:
    for line in text.splitlines():
        cleaned = line.strip()
        if cleaned:
            return re.sub(r"\s+", " ", cleaned)
    return "unknown"


def probe_tool(probe: ToolProbe) -> dict[str, Any]:
    """Resolve one external binary and return provenance metadata."""

    binary = probe.command[0]
    resolved = shutil.which(binary)
    if resolved is None:
        if probe.required:
            raise FileNotFoundError(
                f"required external tool not found on PATH: {binary}",
            )
        return {"available": False, "required": False}

    completed = subprocess.run(
        list(probe.command),
        check=False,
        capture_output=True,
        text=True,
    )
    output = "\n".join(
        part for part in (completed.stdout, completed.stderr) if part.strip()
    )
    if completed.returncode != 0 and probe.required:
        raise RuntimeError(
            f"required external tool failed version probe: "
            f"{' '.join(probe.command)}",
        )

    path = Path(resolved).resolve()
    return {
        "available": completed.returncode == 0,
        "path": str(path),
        "required": probe.required,
        "sha256": file_sha256(path),
        "version": _first_version_line(output),
    }


def resolve_tools(
    probes: tuple[ToolProbe, ...] = DEFAULT_TOOL_PROBES,
) -> dict[str, dict[str, Any]]:
    """Resolve all configured external tools."""

    return {probe.name: probe_tool(probe) for probe in probes}


def write_tools_resolved(
    tools: dict[str, dict[str, Any]],
    out: Path,
) -> Path:
    """Write resolved external-tool provenance as YAML."""

    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(yaml.safe_dump(tools, sort_keys=True), encoding="utf-8")
    return out


def read_tools_resolved(path: Path) -> dict[str, Any]:
    """Read a tools_resolved.yaml file."""

    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"tools resolved file must contain a mapping: {path}")
    return payload
