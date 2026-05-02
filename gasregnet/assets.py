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

Downloader = Literal["auto", "aria2", "wget", "urllib"]
CombineMode = Literal["single", "concat_text"]


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
        combine = raw.get("combine", "single")
        if combine not in {"single", "concat_text"}:
            raise ValueError(f"asset {name} combine must be single or concat_text")
        if combine == "single" and len(urls) != 1:
            raise ValueError(
                f"asset {name} uses single combine but has {len(urls)} URLs",
            )
        sha256 = raw.get("sha256")
        if sha256 is not None and not isinstance(sha256, str):
            raise ValueError(f"asset {name} sha256 must be a string when present")
        assets.append(raw)
    return assets


def _download_with_aria2(url: str, output: Path) -> None:
    subprocess.run(
        [
            "aria2c",
            "--quiet=true",
            "--allow-overwrite=true",
            "--auto-file-renaming=false",
            "--continue=true",
            "--max-connection-per-server=8",
            "--split=8",
            "--min-split-size=1M",
            "--dir",
            str(output.parent),
            "--out",
            output.name,
            url,
        ],
        check=True,
    )


def _download_with_wget(url: str, output: Path) -> None:
    subprocess.run(["wget", "-q", "-O", str(output), url], check=True)


def _download_with_urllib(url: str, output: Path) -> None:
    with urlopen(url, timeout=60) as response:
        output.write_bytes(cast(bytes, response.read()))


def _download(url: str, output: Path, downloader: Downloader) -> None:
    if downloader == "aria2" or (
        downloader == "auto" and shutil.which("aria2c") is not None
    ):
        _download_with_aria2(url, output)
        return
    if downloader == "wget" or (
        downloader == "auto" and shutil.which("wget") is not None
    ):
        _download_with_wget(url, output)
        return
    _download_with_urllib(url, output)


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _combine_downloads(paths: list[Path], combine: CombineMode) -> bytes:
    if combine == "single":
        return paths[0].read_bytes()
    return b"".join(
        payload if payload.endswith(b"\n") else payload + b"\n"
        for payload in (path.read_bytes() for path in paths)
    )


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
            raise ValueError(
                f"asset {asset['name']} hash mismatch: "
                f"expected {expected_hash}, observed {existing_hash}",
            )

        combine = cast(CombineMode, asset.get("combine", "single"))
        with tempfile.TemporaryDirectory(prefix=f"gasregnet-{asset['name']}-") as tmp:
            temp_dir = Path(tmp)
            downloads = [
                temp_dir / f"part-{index}" for index, _ in enumerate(asset["urls"])
            ]
            for url, download in zip(asset["urls"], downloads, strict=True):
                _download(url, download, downloader)
            payload = _combine_downloads(downloads, combine)
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
