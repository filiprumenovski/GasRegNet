"""Fetch versioned external assets declared in configs/assets.yaml."""

from __future__ import annotations

import argparse
from pathlib import Path

from gasregnet.assets import fetch_assets


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("configs/assets.yaml"),
        help="YAML asset manifest.",
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=Path("."),
        help="Repository root for relative output paths.",
    )
    parser.add_argument(
        "--downloader",
        choices=["auto", "aria2", "wget", "urllib"],
        default="auto",
        help="Download backend. auto prefers aria2c, then wget.",
    )
    parser.add_argument("--force", action="store_true", help="Refetch existing assets.")
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    written = fetch_assets(
        args.manifest,
        root=args.root,
        downloader=args.downloader,
        force=args.force,
    )
    for path in written:
        print(path)


if __name__ == "__main__":
    main()
