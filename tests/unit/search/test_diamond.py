from __future__ import annotations

from pathlib import Path

import pytest

from gasregnet.errors import MissingInputError
from gasregnet.search.diamond import parse_diamond_output, run_diamond


def test_parse_diamond_output_reads_typed_columns(tmp_path: Path) -> None:
    output = tmp_path / "hits.tsv"
    output.write_text(
        "q1\ts1\t95.5\t100\t1\t0\t1\t100\t5\t104\t1e-20\t80.0\t99.0\t98.0\n",
        encoding="utf-8",
    )

    hits = parse_diamond_output(output)

    assert hits.height == 1
    assert hits["query_id"].item() == "q1"
    assert hits["length"].item() == 100
    assert hits["evalue"].item() == pytest.approx(1e-20)


def test_parse_diamond_output_rejects_missing_file(tmp_path: Path) -> None:
    with pytest.raises(MissingInputError):
        parse_diamond_output(tmp_path / "missing.tsv")


def test_run_diamond_rejects_missing_query_before_binary_lookup(tmp_path: Path) -> None:
    with pytest.raises(MissingInputError, match="query FASTA"):
        run_diamond(
            tmp_path / "missing.faa",
            tmp_path / "missing.dmnd",
            tmp_path / "hits.tsv",
        )


def test_run_diamond_rejects_missing_binary(tmp_path: Path, monkeypatch) -> None:  # type: ignore[no-untyped-def]
    query = tmp_path / "query.faa"
    db = tmp_path / "db.dmnd"
    query.write_text(">q\nM\n", encoding="utf-8")
    db.write_text("", encoding="utf-8")
    monkeypatch.setenv("PATH", "")

    with pytest.raises(MissingInputError, match="diamond binary"):
        run_diamond(query, db, tmp_path / "hits.tsv")
