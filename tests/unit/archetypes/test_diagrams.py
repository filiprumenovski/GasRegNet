from __future__ import annotations

from pathlib import Path

from gasregnet.archetypes.cluster import cluster_archetypes
from gasregnet.archetypes.diagrams import render_archetype_diagrams
from tests.unit.archetypes.test_cluster import _candidates, _loci


def test_render_archetype_diagrams_writes_svg_and_png(tmp_path: Path) -> None:
    archetypes = cluster_archetypes(_loci(), _candidates())

    outputs = render_archetype_diagrams(archetypes, tmp_path)

    assert len(outputs) == archetypes.height
    for svg_path, png_path in outputs:
        assert svg_path.exists()
        assert png_path.exists()
        assert svg_path.read_text(encoding="utf-8").lstrip().startswith("<?xml")
        assert png_path.stat().st_size > 0
