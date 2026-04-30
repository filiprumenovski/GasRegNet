"""Multiple-sequence alignment helpers."""

from __future__ import annotations

from pathlib import Path

import polars as pl

from gasregnet.io.fasta import read_fasta


def read_alignment_fasta(path: Path) -> pl.DataFrame:
    """Read an aligned FASTA into sequence rows."""

    records = list(read_fasta(path))
    return pl.DataFrame(
        {
            "sequence_id": [record[0] for record in records],
            "aligned_sequence": [record[2] for record in records],
        },
        schema_overrides={"sequence_id": pl.Utf8, "aligned_sequence": pl.Utf8},
    )


def conserved_residues(
    alignment: pl.DataFrame,
    *,
    min_fraction: float = 1.0,
) -> pl.DataFrame:
    """Return alignment columns conserved at or above ``min_fraction``."""

    if alignment.is_empty():
        return pl.DataFrame(
            schema={
                "alignment_position": pl.Int64,
                "residue": pl.Utf8,
                "conservation_fraction": pl.Float64,
            },
        )
    sequences = [str(value) for value in alignment["aligned_sequence"].to_list()]
    width = max(len(sequence) for sequence in sequences)
    rows: list[dict[str, object]] = []
    for index in range(width):
        residues = [
            sequence[index]
            for sequence in sequences
            if index < len(sequence) and sequence[index] != "-"
        ]
        if not residues:
            continue
        residue = sorted(
            set(residues),
            key=lambda value: (-residues.count(value), value),
        )[0]
        fraction = residues.count(residue) / len(residues)
        if fraction >= min_fraction:
            rows.append(
                {
                    "alignment_position": index + 1,
                    "residue": residue,
                    "conservation_fraction": fraction,
                },
            )
    return pl.DataFrame(
        rows,
        schema_overrides={
            "alignment_position": pl.Int64,
            "residue": pl.Utf8,
            "conservation_fraction": pl.Float64,
        },
    )
