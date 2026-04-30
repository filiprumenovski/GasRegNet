from __future__ import annotations

from pathlib import Path

import pytest

from gasregnet.errors import ConfigError, MissingInputError
from gasregnet.neighborhoods.retrieve import retrieve_from_efi_sqlite
from scripts.build_test_fixtures import MINI_EFI, build_mini_efi


@pytest.fixture(scope="module", autouse=True)
def mini_efi_fixture() -> None:
    build_mini_efi()


def test_retrieve_from_efi_sqlite_concatenates_analytes() -> None:
    loci, genes = retrieve_from_efi_sqlite(MINI_EFI, ["CO", "CN"])

    assert loci.height == 2
    assert genes.height == 8
    assert set(loci["analyte"].to_list()) == {"CO", "CN"}


def test_retrieve_from_efi_sqlite_applies_cluster_filter() -> None:
    loci, genes = retrieve_from_efi_sqlite(MINI_EFI, ["CO"], cluster_filter=[999])

    assert loci.height == 0
    assert genes.height == 0


def test_retrieve_from_efi_sqlite_requires_analyte() -> None:
    with pytest.raises(ConfigError, match="at least one analyte"):
        retrieve_from_efi_sqlite(MINI_EFI, [])


def test_retrieve_from_efi_sqlite_rejects_missing_file(tmp_path: Path) -> None:
    with pytest.raises(MissingInputError):
        retrieve_from_efi_sqlite(tmp_path / "missing.sqlite", ["CO"])
