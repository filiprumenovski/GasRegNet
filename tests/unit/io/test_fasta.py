from __future__ import annotations

from pathlib import Path

import pytest

from gasregnet.errors import MissingInputError, SchemaError
from gasregnet.io.fasta import read_fasta


def test_read_fasta_yields_accession_description_sequence(tmp_path: Path) -> None:
    path = tmp_path / "seqs.faa"
    path.write_text(">seq1 description here\nMAGA\n>seq2\nTT\n", encoding="utf-8")

    records = list(read_fasta(path))

    assert records == [
        ("seq1", "seq1 description here", "MAGA"),
        ("seq2", "seq2", "TT"),
    ]


def test_read_fasta_rejects_missing_file(tmp_path: Path) -> None:
    with pytest.raises(MissingInputError):
        list(read_fasta(tmp_path / "missing.faa"))


def test_read_fasta_rejects_empty_file(tmp_path: Path) -> None:
    path = tmp_path / "empty.faa"
    path.write_text("", encoding="utf-8")

    with pytest.raises(SchemaError, match="no records"):
        list(read_fasta(path))
