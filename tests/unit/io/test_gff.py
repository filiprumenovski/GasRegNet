from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from gasregnet.errors import SchemaError
from gasregnet.io.gff import parse_attributes, read_gff3


def test_parse_attributes_decodes_gff3_values() -> None:
    assert parse_attributes("ID=gene1;Name=coxL%20gene") == {
        "ID": "gene1",
        "Name": "coxL gene",
    }


def test_read_gff3_accepts_version_3(tmp_path: Path) -> None:
    path = tmp_path / "genes.gff3"
    path.write_text(
        "##gff-version 3\n"
        "contig\tRefSeq\tCDS\t1\t99\t.\t+\t0\tID=cds1;Name=coxL\n",
        encoding="utf-8",
    )

    frame = read_gff3(path)

    assert frame.height == 1
    assert frame["seqid"].item() == "contig"
    assert frame["attributes"].item()["ID"] == "cds1"


def test_read_gff3_rejects_gff2(tmp_path: Path) -> None:
    path = tmp_path / "genes.gff"
    path.write_text("##gff-version 2\n", encoding="utf-8")

    with pytest.raises(SchemaError, match="GFF3"):
        read_gff3(path)


def test_read_gff3_accepts_gzip(tmp_path: Path) -> None:
    path = tmp_path / "genes.gff3.gz"
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(
            "##gff-version 3\n"
            "contig\tRefSeq\tCDS\t1\t99\t.\t+\t0\tID=cds1;Name=coxL\n",
        )

    frame = read_gff3(path)

    assert frame.height == 1
    assert frame["feature_type"].item() == "CDS"
