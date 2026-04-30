from __future__ import annotations

from pathlib import Path

import pytest

from gasregnet.errors import MissingInputError, SchemaError
from gasregnet.io.parquet import read_parquet, write_parquet
from gasregnet.schemas import LociSchema
from tests.unit.test_schemas import loci_frame


def test_parquet_round_trip_preserves_schema(tmp_path: Path) -> None:
    path = tmp_path / "loci.parquet"

    write_parquet(loci_frame(), path, LociSchema)
    loaded = read_parquet(path, LociSchema)

    assert loaded.equals(loci_frame())


def test_read_parquet_rejects_missing_file(tmp_path: Path) -> None:
    with pytest.raises(MissingInputError):
        read_parquet(tmp_path / "missing.parquet", LociSchema)


def test_write_parquet_validates_before_writing(tmp_path: Path) -> None:
    with pytest.raises(SchemaError):
        write_parquet(
            loci_frame().drop("locus_id"),
            tmp_path / "bad.parquet",
            LociSchema,
        )
