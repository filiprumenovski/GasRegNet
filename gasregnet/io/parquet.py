"""Validated Parquet I/O helpers."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import polars as pl

from gasregnet.errors import MissingInputError
from gasregnet.schemas import DataFrameSchema, validate


def write_parquet(
    df: pl.DataFrame,
    path: Path,
    schema: DataFrameSchema,
    **kwargs: Any,
) -> Path:
    """Validate and write a DataFrame as Parquet."""

    path.parent.mkdir(parents=True, exist_ok=True)
    validated = validate(df, schema)
    validated.write_parquet(path, **kwargs)
    return path


def read_parquet(path: Path, schema: DataFrameSchema) -> pl.DataFrame:
    """Read a Parquet file and validate it against a schema."""

    if not path.exists():
        raise MissingInputError(f"Parquet file does not exist: {path}")
    return validate(pl.read_parquet(path), schema)
