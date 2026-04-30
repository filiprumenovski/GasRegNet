from __future__ import annotations

import hashlib
from pathlib import Path

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
