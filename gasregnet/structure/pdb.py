"""PDB alignment helpers."""

from __future__ import annotations

from pathlib import Path

import polars as pl

from gasregnet.errors import MissingInputError

THREE_TO_ONE = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}


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
    """Map model residues to homolog residues by sequence alignment of CA atoms."""

    model = _ca_residues(model_pdb)
    homolog = _ca_residues(homolog_pdb)
    if chain_id is not None:
        model = [row for row in model if row["chain_id"] == chain_id]
        homolog = [row for row in homolog if row["chain_id"] == chain_id]
    aligned_pairs = _needleman_wunsch_pairs(
        "".join(_one_letter(row) for row in model),
        "".join(_one_letter(row) for row in homolog),
    )
    rows = []
    for index, (model_index, homolog_index) in enumerate(aligned_pairs):
        if model_index is None or homolog_index is None:
            continue
        model_row = model[model_index]
        homolog_row = homolog[homolog_index]
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


def _one_letter(row: dict[str, object]) -> str:
    return THREE_TO_ONE.get(str(row["residue_name"]).upper(), "X")


def _needleman_wunsch_pairs(
    left: str,
    right: str,
    *,
    match: int = 2,
    mismatch: int = -1,
    gap: int = -2,
) -> list[tuple[int | None, int | None]]:
    """Return global sequence-alignment residue index pairs."""

    n = len(left)
    m = len(right)
    scores = [[0] * (m + 1) for _ in range(n + 1)]
    trace = [[""] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        scores[i][0] = scores[i - 1][0] + gap
        trace[i][0] = "up"
    for j in range(1, m + 1):
        scores[0][j] = scores[0][j - 1] + gap
        trace[0][j] = "left"
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diagonal = scores[i - 1][j - 1] + (
                match if left[i - 1] == right[j - 1] else mismatch
            )
            up = scores[i - 1][j] + gap
            left_score = scores[i][j - 1] + gap
            best = max(diagonal, up, left_score)
            scores[i][j] = best
            if best == diagonal:
                trace[i][j] = "diagonal"
            elif best == up:
                trace[i][j] = "up"
            else:
                trace[i][j] = "left"

    pairs: list[tuple[int | None, int | None]] = []
    i = n
    j = m
    while i > 0 or j > 0:
        direction = trace[i][j]
        if direction == "diagonal":
            pairs.append((i - 1, j - 1))
            i -= 1
            j -= 1
        elif direction == "up":
            pairs.append((i - 1, None))
            i -= 1
        else:
            pairs.append((None, j - 1))
            j -= 1
    return list(reversed(pairs))
