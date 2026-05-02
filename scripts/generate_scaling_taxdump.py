"""Generate a minimal deterministic NCBI taxdump for scaling probes."""

from __future__ import annotations

import argparse
from pathlib import Path


def _node(taxon_id: int, parent_id: int, rank: str) -> str:
    return f"{taxon_id}\t|\t{parent_id}\t|\t{rank}\t|\t\t|\n"


def _name(taxon_id: int, name: str) -> str:
    return f"{taxon_id}\t|\t{name}\t|\t\t|\tscientific name\t|\n"


def generate_scaling_taxdump(out_dir: Path, *, datasets: int = 3) -> Path:
    """Write nodes.dmp and names.dmp covering generated scaling taxon IDs."""

    out_dir.mkdir(parents=True, exist_ok=True)
    nodes = [
        _node(1, 1, "no rank"),
        _node(2, 1, "superkingdom"),
        _node(1224, 2, "phylum"),
        _node(1236, 1224, "class"),
        _node(91061, 2, "phylum"),
        _node(1385, 91061, "class"),
    ]
    names = [
        _name(1, "root"),
        _name(2, "Bacteria"),
        _name(1224, "Pseudomonadota"),
        _name(1236, "Alphaproteobacteria"),
        _name(91061, "Bacillota"),
        _name(1385, "Bacilli"),
    ]
    for dataset_index in range(datasets):
        class_id = 1236 if dataset_index % 2 == 0 else 1385
        order_id = 910_000 + dataset_index * 4
        family_id = order_id + 1
        genus_id = order_id + 2
        species_id = 900_000 + dataset_index
        nodes.extend(
            [
                _node(order_id, class_id, "order"),
                _node(family_id, order_id, "family"),
                _node(genus_id, family_id, "genus"),
                _node(species_id, genus_id, "species"),
            ],
        )
        names.extend(
            [
                _name(order_id, f"ScalingOrder{dataset_index}"),
                _name(family_id, f"ScalingFamily{dataset_index}"),
                _name(genus_id, f"ScalingGenus{dataset_index}"),
                _name(species_id, f"Scaling species {dataset_index}"),
            ],
        )

    (out_dir / "nodes.dmp").write_text("".join(nodes), encoding="utf-8")
    (out_dir / "names.dmp").write_text("".join(names), encoding="utf-8")
    return out_dir


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--datasets", type=int, default=3)
    args = parser.parse_args()
    print(generate_scaling_taxdump(args.out_dir, datasets=args.datasets))


if __name__ == "__main__":
    main()
