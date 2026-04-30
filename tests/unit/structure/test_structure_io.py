from __future__ import annotations

import json
from pathlib import Path

import pytest

from gasregnet.errors import MissingInputError
from gasregnet.structure.alphafold import read_pae_json, read_plddt_from_pdb
from gasregnet.structure.msa import conserved_residues, read_alignment_fasta
from gasregnet.structure.pdb import residue_mapping_by_order


def _atom_line(
    serial: int,
    residue: str,
    residue_number: int,
    *,
    chain: str = "A",
    bfactor: float = 90.0,
) -> str:
    return (
        f"ATOM  {serial:5d}  CA  {residue:>3s} {chain}{residue_number:4d}    "
        f"{float(serial):8.3f}{0.0:8.3f}{0.0:8.3f}{1.0:6.2f}{bfactor:6.2f}"
        "           C\n"
    )


def test_read_alignment_fasta_and_conserved_residues(tmp_path: Path) -> None:
    fasta = tmp_path / "alignment.faa"
    fasta.write_text(">a\nMA-\n>b\nMAA\n", encoding="utf-8")

    alignment = read_alignment_fasta(fasta)
    conserved = conserved_residues(alignment, min_fraction=1.0)

    assert alignment.height == 2
    assert conserved["alignment_position"].to_list() == [1, 2, 3]
    assert conserved["residue"].to_list() == ["M", "A", "A"]


def test_read_plddt_from_pdb_uses_residue_bfactor(tmp_path: Path) -> None:
    pdb = tmp_path / "model.pdb"
    pdb.write_text(
        _atom_line(1, "MET", 1, bfactor=91.5)
        + _atom_line(2, "ALA", 2, bfactor=70.0),
        encoding="utf-8",
    )

    confidence = read_plddt_from_pdb(pdb)

    assert confidence.height == 2
    assert confidence["plddt"].to_list() == [91.5, 70.0]


def test_read_pae_json_flattens_error_matrix(tmp_path: Path) -> None:
    pae = tmp_path / "pae.json"
    pae.write_text(
        json.dumps({"predicted_aligned_error": [[0.0, 1.5], [1.5, 0.0]]}),
        encoding="utf-8",
    )

    frame = read_pae_json(pae)

    assert frame.height == 4
    assert frame.filter(frame["residue_i"] == 1)["pae"].to_list() == [0.0, 1.5]


def test_residue_mapping_by_order_pairs_ca_atoms(tmp_path: Path) -> None:
    model = tmp_path / "model.pdb"
    homolog = tmp_path / "homolog.pdb"
    model.write_text(
        _atom_line(1, "MET", 10) + _atom_line(2, "ALA", 11),
        encoding="utf-8",
    )
    homolog.write_text(
        _atom_line(1, "MET", 50) + _atom_line(2, "SER", 51),
        encoding="utf-8",
    )

    mapping = residue_mapping_by_order(model, homolog)

    assert mapping["model_residue_number"].to_list() == [10, 11]
    assert mapping["homolog_residue_number"].to_list() == [50, 51]


def test_structure_parsers_reject_missing_files(tmp_path: Path) -> None:
    with pytest.raises(MissingInputError):
        read_plddt_from_pdb(tmp_path / "missing.pdb")
