"""FASTA parsing helpers."""

from __future__ import annotations

from collections.abc import Iterator
from pathlib import Path

from Bio import SeqIO

from gasregnet.errors import MissingInputError, SchemaError


def read_fasta(path: Path) -> Iterator[tuple[str, str, str]]:
    """Yield ``(accession, description, sequence)`` tuples from a FASTA file."""

    if not path.exists():
        raise MissingInputError(f"FASTA file does not exist: {path}")

    n_records = 0
    for record in SeqIO.parse(path, "fasta"):  # type: ignore[no-untyped-call]
        n_records += 1
        sequence = str(record.seq)
        if not sequence:
            raise SchemaError(f"FASTA record {record.id!r} has an empty sequence")
        yield record.id, record.description, sequence

    if n_records == 0:
        raise SchemaError(f"FASTA file contains no records: {path}")
