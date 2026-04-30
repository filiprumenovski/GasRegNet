"""FASTA parsing helpers."""

from __future__ import annotations

import gzip
from collections.abc import Iterator
from pathlib import Path
from typing import TextIO

from Bio import SeqIO

from gasregnet.errors import MissingInputError, SchemaError


def _open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def read_fasta(path: Path) -> Iterator[tuple[str, str, str]]:
    """Yield ``(accession, description, sequence)`` tuples from a FASTA file."""

    if not path.exists():
        raise MissingInputError(f"FASTA file does not exist: {path}")

    n_records = 0
    with _open_text(path) as handle:
        for record in SeqIO.parse(handle, "fasta"):  # type: ignore[no-untyped-call]
            n_records += 1
            sequence = str(record.seq)
            if not sequence:
                raise SchemaError(f"FASTA record {record.id!r} has an empty sequence")
            yield record.id, record.description, sequence

    if n_records == 0:
        raise SchemaError(f"FASTA file contains no records: {path}")
