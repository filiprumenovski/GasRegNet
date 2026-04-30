"""PDB alignment helpers."""

from __future__ import annotations

from pathlib import Path

import polars as pl

from gasregnet.errors import MissingInputError


def _ca_residues(path: Path) -> list[dict[str, object]]:
    if not path.exists():
        raise MissingInputError(f"PDB file does not exist: {path}")
    residues: list[dict[str, object]] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("ATOM") or line[12:16].strip() != "CA":
                continue
            residues.append(
                {
                    "chain_id": line[21].strip() or "_",
                    "residue_number": int(line[22:26]),
                    "insertion_code": line[26].strip(),
                    "residue_name": line[17:20].strip(),
                    "x": float(line[30:38]),
                    "y": float(line[38:46]),
                    "z": float(line[46:54]),
                },
            )
    return residues


def residue_mapping_by_order(
    model_pdb: Path,
    homolog_pdb: Path,
    *,
    chain_id: str | None = None,
) -> pl.DataFrame:
    """Map model residues to homolog residues by ordered CA atoms."""

    model = _ca_residues(model_pdb)
    homolog = _ca_residues(homolog_pdb)
    if chain_id is not None:
        model = [row for row in model if row["chain_id"] == chain_id]
        homolog = [row for row in homolog if row["chain_id"] == chain_id]
    rows = []
    for index, (model_row, homolog_row) in enumerate(zip(model, homolog, strict=False)):
        rows.append(
            {
                "alignment_index": index + 1,
                "model_chain_id": model_row["chain_id"],
                "model_residue_number": model_row["residue_number"],
                "model_residue_name": model_row["residue_name"],
                "homolog_chain_id": homolog_row["chain_id"],
                "homolog_residue_number": homolog_row["residue_number"],
                "homolog_residue_name": homolog_row["residue_name"],
            },
        )
    return pl.DataFrame(
        rows,
        schema_overrides={
            "alignment_index": pl.Int64,
            "model_chain_id": pl.Utf8,
            "model_residue_number": pl.Int64,
            "model_residue_name": pl.Utf8,
            "homolog_chain_id": pl.Utf8,
            "homolog_residue_number": pl.Int64,
            "homolog_residue_name": pl.Utf8,
        },
    )
