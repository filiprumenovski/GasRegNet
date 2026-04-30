from __future__ import annotations

from pathlib import Path

import pytest

from gasregnet.errors import MissingInputError, SchemaError
from gasregnet.io.sqlite_efi import read_efi_sqlite
from scripts.build_test_fixtures import MINI_EFI, build_mini_efi


@pytest.fixture(scope="module", autouse=True)
def mini_efi_fixture() -> None:
    build_mini_efi()


def test_read_efi_sqlite_returns_valid_loci_and_genes() -> None:
    loci, genes = read_efi_sqlite(MINI_EFI, "CO")

    assert loci.height == 1
    assert genes.height == 4
    assert loci["locus_id"].item() == "CO_7_COX_L_ANCHOR"


def test_read_efi_sqlite_computes_relative_indices() -> None:
    _, genes = read_efi_sqlite(MINI_EFI, "CO")

    relative = dict(zip(genes["gene_accession"], genes["relative_index"], strict=True))
    assert relative == {
        "REG_UP": -2,
        "COX_M": -1,
        "COX_L_ANCHOR": 0,
        "COX_S": 1,
    }
    anchor = genes.filter(genes["gene_accession"] == "COX_L_ANCHOR")
    assert anchor["is_anchor"].item() is True
    assert anchor["relative_start"].item() == 0


def test_read_efi_sqlite_applies_cluster_filter() -> None:
    loci, genes = read_efi_sqlite(MINI_EFI, "CO", cluster_filter=[999])

    assert loci.height == 0
    assert genes.height == 0


def test_read_efi_sqlite_rejects_missing_file(tmp_path: Path) -> None:
    with pytest.raises(MissingInputError):
        read_efi_sqlite(tmp_path / "missing.sqlite", "CO")


def test_read_efi_sqlite_rejects_missing_tables(tmp_path: Path) -> None:
    bad_sqlite = tmp_path / "bad.sqlite"
    bad_sqlite.write_bytes(b"")

    with pytest.raises(SchemaError, match="missing required table"):
        read_efi_sqlite(bad_sqlite, "CO")
