"""External asset manifest fetching."""

from __future__ import annotations

import hashlib
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any, Literal, cast
from urllib.request import urlopen

import yaml  # type: ignore[import-untyped]

Downloader = Literal["auto", "wget", "urllib"]


def _read_manifest(path: Path) -> list[dict[str, Any]]:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict) or not isinstance(payload.get("assets"), list):
        raise ValueError(f"asset manifest must contain an assets list: {path}")
    assets: list[dict[str, Any]] = []
    for raw in payload["assets"]:
        if not isinstance(raw, dict):
            raise ValueError("asset entries must be mappings")
        name = raw.get("name")
        output = raw.get("output")
        urls = raw.get("urls")
        if not isinstance(name, str) or not name:
            raise ValueError("asset entry is missing string field: name")
        if not isinstance(output, str) or not output:
            raise ValueError(f"asset {name} is missing string field: output")
        if not isinstance(urls, list) or not all(isinstance(url, str) for url in urls):
            raise ValueError(f"asset {name} is missing string list field: urls")
        sha256 = raw.get("sha256")
        if sha256 is not None and not isinstance(sha256, str):
            raise ValueError(f"asset {name} sha256 must be a string when present")
        assets.append(raw)
    return assets


def _download_with_wget(url: str) -> bytes:
    result = subprocess.run(
        ["wget", "-q", "-O", "-", url],
        check=True,
        stdout=subprocess.PIPE,
    )
    return result.stdout


def _download_with_urllib(url: str) -> bytes:
    with urlopen(url, timeout=60) as response:
        return cast(bytes, response.read())


def _download(url: str, downloader: Downloader) -> bytes:
    if downloader == "wget" or (
        downloader == "auto" and shutil.which("wget") is not None
    ):
        return _download_with_wget(url)
    return _download_with_urllib(url)


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def fetch_assets(
    manifest: Path,
    *,
    root: Path,
    downloader: Downloader = "auto",
    force: bool = False,
) -> list[Path]:
    """Fetch manifest assets and return written output paths."""

    written: list[Path] = []
    for asset in _read_manifest(manifest):
        output = root / str(asset["output"])
        if output.exists() and not force:
            existing_hash = _sha256(output.read_bytes())
            expected_hash = asset.get("sha256")
            if expected_hash is None or existing_hash == expected_hash:
                written.append(output)
                continue

        chunks = [_download(url, downloader) for url in asset["urls"]]
        payload = b"".join(
            chunk if chunk.endswith(b"\n") else chunk + b"\n" for chunk in chunks
        )
        expected_hash = asset.get("sha256")
        observed_hash = _sha256(payload)
        if expected_hash is not None and observed_hash != expected_hash:
            raise ValueError(
                f"asset {asset['name']} hash mismatch: "
                f"expected {expected_hash}, observed {observed_hash}",
            )

        output.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(
            dir=output.parent,
            prefix=f".{output.name}.",
            delete=False,
        ) as handle:
            handle.write(payload)
            temp_path = Path(handle.name)
        temp_path.replace(output)
        written.append(output)
    return written
