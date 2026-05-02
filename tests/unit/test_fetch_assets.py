from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

from gasregnet.assets import fetch_assets


def test_fetch_assets_concatenates_urls_and_checks_hash(tmp_path: Path) -> None:
    source_a = tmp_path / "a.faa"
    source_b = tmp_path / "b.faa"
    source_a.write_text(">a\nMA\n", encoding="utf-8")
    source_b.write_text(">b\nMG", encoding="utf-8")
    payload = b">a\nMA\n>b\nMG\n"
    digest = hashlib.sha256(payload).hexdigest()
    manifest = tmp_path / "assets.yaml"
    manifest.write_text(
        f"""
assets:
  - name: mini
    output: data/seeds/mini.faa
    combine: concat_text
    sha256: {digest}
    urls:
      - {source_a.as_uri()}
      - {source_b.as_uri()}
""",
        encoding="utf-8",
    )

    written = fetch_assets(
        manifest,
        root=tmp_path,
        downloader="urllib",
        force=True,
    )

    assert written == [tmp_path / "data" / "seeds" / "mini.faa"]
    assert written[0].read_bytes() == payload


def test_fetch_assets_preserves_single_file_bytes(tmp_path: Path) -> None:
    source = tmp_path / "archive.gz"
    payload = b"\x1f\x8b\x08\x00binary"
    source.write_bytes(payload)
    digest = hashlib.sha256(payload).hexdigest()
    manifest = tmp_path / "assets.yaml"
    manifest.write_text(
        f"""
assets:
  - name: binary
    output: data/external/archive.gz
    sha256: {digest}
    urls:
      - {source.as_uri()}
""",
        encoding="utf-8",
    )

    written = fetch_assets(
        manifest,
        root=tmp_path,
        downloader="urllib",
        force=True,
    )

    assert written == [tmp_path / "data" / "external" / "archive.gz"]
    assert written[0].read_bytes() == payload


def test_fetch_assets_fails_on_sha256_mismatch(tmp_path: Path) -> None:
    source = tmp_path / "seed.faa"
    source.write_text(">seed\nMA\n", encoding="utf-8")
    manifest = tmp_path / "assets.yaml"
    manifest.write_text(
        f"""
assets:
  - name: bad-seed
    output: data/seeds/bad.faa
    sha256: '{"0" * 64}'
    urls:
      - {source.as_uri()}
""",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="bad-seed hash mismatch"):
        fetch_assets(
            manifest,
            root=tmp_path,
            downloader="urllib",
            force=True,
        )


def test_fetch_assets_fails_on_existing_file_sha256_mismatch(tmp_path: Path) -> None:
    source = tmp_path / "source.faa"
    source.write_text(">seed\nMA\n", encoding="utf-8")
    expected = hashlib.sha256(source.read_bytes()).hexdigest()
    output = tmp_path / "data" / "seeds" / "seed.faa"
    output.parent.mkdir(parents=True)
    output.write_text(">corrupt\nXX\n", encoding="utf-8")
    manifest = tmp_path / "assets.yaml"
    manifest.write_text(
        f"""
assets:
  - name: cached-seed
    output: data/seeds/seed.faa
    sha256: {expected}
    urls:
      - {source.as_uri()}
""",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="cached-seed hash mismatch"):
        fetch_assets(
            manifest,
            root=tmp_path,
            downloader="urllib",
            force=False,
        )
