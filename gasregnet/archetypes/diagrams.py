"""Archetype diagram rendering."""

from __future__ import annotations

import re
from pathlib import Path

import matplotlib.pyplot as plt
import polars as pl
from matplotlib.patches import FancyArrow

TOKEN_RE = re.compile(r"\[(-?\d+):([^\]]+)\]")


def _tokens(architecture: str) -> list[tuple[int, str]]:
    return [
        (int(match.group(1)), match.group(2))
        for match in TOKEN_RE.finditer(architecture)
    ]


def _safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value)


def render_archetype_diagram(
    archetype: dict[str, object],
    out_dir: Path,
) -> tuple[Path, Path]:
    """Render one archetype as SVG and PNG gene-arrow diagrams."""

    out_dir.mkdir(parents=True, exist_ok=True)
    archetype_id = str(archetype["archetype_id"])
    architecture = str(archetype["architecture_string"])
    tokens = _tokens(architecture)
    if not tokens:
        tokens = [(0, "empty")]

    fig_width = max(4.0, 1.2 * len(tokens))
    fig, ax = plt.subplots(figsize=(fig_width, 1.6))
    ax.set_axis_off()
    ax.set_xlim(-0.8, len(tokens) - 0.2)
    ax.set_ylim(-0.8, 0.8)

    for index, (relative_index, label) in enumerate(tokens):
        color = "#4C78A8" if relative_index == 0 else "#F58518"
        ax.add_patch(
            FancyArrow(
                index - 0.35,
                0,
                0.7,
                0,
                width=0.22,
                head_width=0.38,
                head_length=0.16,
                length_includes_head=True,
                color=color,
            ),
        )
        ax.text(index, 0.42, str(relative_index), ha="center", va="bottom", fontsize=8)
        ax.text(index, -0.42, label, ha="center", va="top", fontsize=7)

    ax.set_title(
        f"{archetype_id}: n={archetype['n_loci']} loci",
        fontsize=9,
        loc="left",
    )
    fig.tight_layout()

    base = out_dir / _safe_name(archetype_id)
    svg_path = base.with_suffix(".svg")
    png_path = base.with_suffix(".png")
    fig.savefig(svg_path)
    fig.savefig(png_path, dpi=300)
    plt.close(fig)
    return svg_path, png_path


def render_archetype_diagrams(
    archetypes: pl.DataFrame,
    out_dir: Path,
) -> list[tuple[Path, Path]]:
    """Render SVG and PNG diagrams for every archetype."""

    return [
        render_archetype_diagram(archetype, out_dir)
        for archetype in archetypes.iter_rows(named=True)
    ]
