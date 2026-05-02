from __future__ import annotations

import stat
from pathlib import Path

import pytest

from gasregnet.check_tools import ToolProbe, resolve_tools, write_tools_resolved
from gasregnet.hashing import file_sha256


def _write_executable(path: Path, body: str) -> Path:
    path.write_text(body, encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IXUSR)
    return path


def test_resolve_tools_records_path_version_and_binary_hash(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    binary = _write_executable(
        tmp_path / "fake-tool",
        "#!/bin/sh\necho fake-tool 1.2.3\n",
    )
    monkeypatch.setenv("PATH", str(tmp_path))

    resolved = resolve_tools((ToolProbe("fake", ("fake-tool", "--version")),))

    assert resolved["fake"]["available"] is True
    assert resolved["fake"]["path"] == str(binary.resolve())
    assert resolved["fake"]["sha256"] == file_sha256(binary.resolve())
    assert resolved["fake"]["version"] == "fake-tool 1.2.3"


def test_required_missing_tool_fails(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("PATH", "")

    with pytest.raises(FileNotFoundError, match="required external tool"):
        resolve_tools((ToolProbe("missing", ("definitely-missing",)),))


def test_optional_missing_tool_is_recorded(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("PATH", "")

    resolved = resolve_tools(
        (ToolProbe("optional", ("definitely-missing",), required=False),),
    )

    assert resolved == {"optional": {"available": False, "required": False}}


def test_write_tools_resolved_outputs_yaml(tmp_path: Path) -> None:
    out = write_tools_resolved(
        {"fake": {"available": True, "version": "fake 1"}},
        tmp_path / "tools_resolved.yaml",
    )

    assert out.read_text(encoding="utf-8") == (
        "fake:\n  available: true\n  version: fake 1\n"
    )
