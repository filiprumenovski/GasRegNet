"""Path helpers."""

from __future__ import annotations

from pathlib import Path


def ensure_out_dir(out_dir: Path) -> Path:
    """Create a standard run output directory and return it."""

    out_dir.mkdir(parents=True, exist_ok=True)
    for child in ("logs", "tables", "figures", "intermediate"):
        (out_dir / child).mkdir(exist_ok=True)
    return out_dir
