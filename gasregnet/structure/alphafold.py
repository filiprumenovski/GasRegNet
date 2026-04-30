"""AlphaFold output ingestion."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import polars as pl

from gasregnet.errors import MissingInputError, SchemaError


def read_plddt_from_pdb(pdb_path: Path, *, chain_id: str | None = None) -> pl.DataFrame:
    """Read per-residue pLDDT values from AlphaFold PDB B-factors."""

    if not pdb_path.exists():
        raise MissingInputError(f"AlphaFold PDB does not exist: {pdb_path}")

    rows: dict[tuple[str, int, str], dict[str, object]] = {}
    with pdb_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue
            chain = line[21].strip() or "_"
            if chain_id is not None and chain != chain_id:
                continue
            residue_number = int(line[22:26])
            insertion_code = line[26].strip()
            key = (chain, residue_number, insertion_code)
            rows.setdefault(
                key,
                {
                    "chain_id": chain,
                    "residue_number": residue_number,
                    "insertion_code": insertion_code,
                    "residue_name": line[17:20].strip(),
                    "plddt": float(line[60:66]),
                },
            )

    return pl.DataFrame(
        list(rows.values()),
        schema_overrides={
            "chain_id": pl.Utf8,
            "residue_number": pl.Int64,
            "insertion_code": pl.Utf8,
            "residue_name": pl.Utf8,
            "plddt": pl.Float64,
        },
    )


def read_pae_json(path: Path) -> pl.DataFrame:
    """Read AlphaFold predicted aligned error JSON into a long table."""

    if not path.exists():
        raise MissingInputError(f"AlphaFold PAE JSON does not exist: {path}")
    data = json.loads(path.read_text(encoding="utf-8"))
    matrix = _pae_matrix(data)
    rows = [
        {"residue_i": i + 1, "residue_j": j + 1, "pae": float(value)}
        for i, row in enumerate(matrix)
        for j, value in enumerate(row)
    ]
    return pl.DataFrame(
        rows,
        schema_overrides={
            "residue_i": pl.Int64,
            "residue_j": pl.Int64,
            "pae": pl.Float64,
        },
    )


def _pae_matrix(data: Any) -> list[list[float]]:
    if isinstance(data, dict):
        matrix = data.get("predicted_aligned_error") or data.get("pae")
    elif isinstance(data, list) and data and isinstance(data[0], dict):
        matrix = data[0].get("predicted_aligned_error") or data[0].get("pae")
    else:
        matrix = data
    if not isinstance(matrix, list) or not all(isinstance(row, list) for row in matrix):
        raise SchemaError("PAE JSON must contain a rectangular error matrix")
    return [[float(value) for value in row] for row in matrix]
