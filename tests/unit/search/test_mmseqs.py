from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest

from gasregnet.errors import MissingInputError
from gasregnet.search.mmseqs import cluster_sequences, parse_mmseqs_cluster_tsv


def test_parse_mmseqs_cluster_tsv_marks_representatives(tmp_path: Path) -> None:
    output = tmp_path / "clusters_cluster.tsv"
    output.write_text(
        "rep1\trep1\nrep1\tmember1\nrep2\trep2\n",
        encoding="utf-8",
    )

    clusters = parse_mmseqs_cluster_tsv(output)

    assert clusters.height == 3
    assert clusters.columns == ["sequence_id", "cluster_id", "is_representative"]
    assert clusters.filter(clusters["sequence_id"] == "rep1")[
        "is_representative"
    ].item() is True
    assert clusters.filter(clusters["sequence_id"] == "member1")[
        "cluster_id"
    ].item() == "rep1"


def test_parse_mmseqs_cluster_tsv_rejects_missing_file(tmp_path: Path) -> None:
    with pytest.raises(MissingInputError):
        parse_mmseqs_cluster_tsv(tmp_path / "missing.tsv")


def test_cluster_sequences_rejects_missing_input_before_binary_lookup(
    tmp_path: Path,
) -> None:
    with pytest.raises(MissingInputError, match="input FASTA"):
        cluster_sequences(tmp_path / "missing.faa", tmp_path / "clusters")


def test_cluster_sequences_runs_mmseqs_and_parses_output(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    input_faa = tmp_path / "input.faa"
    input_faa.write_text(">rep1\nM\n>member1\nM\n", encoding="utf-8")

    def fake_run(command: list[str], **_: Any) -> None:
        assert command[1] == "easy-cluster"
        assert command[2] == str(input_faa)
        Path(f"{command[3]}_cluster.tsv").write_text(
            "rep1\trep1\nrep1\tmember1\n",
            encoding="utf-8",
        )

    monkeypatch.setenv("PATH", "/usr/bin")
    monkeypatch.setattr("shutil.which", lambda _: "/usr/bin/mmseqs")
    monkeypatch.setattr("subprocess.run", fake_run)

    clusters = cluster_sequences(input_faa, tmp_path / "clusters")

    assert clusters.height == 2
    assert clusters["cluster_id"].to_list() == ["rep1", "rep1"]


def test_cluster_sequences_rejects_missing_binary(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    input_faa = tmp_path / "input.faa"
    input_faa.write_text(">seq\nM\n", encoding="utf-8")
    monkeypatch.setenv("PATH", "")

    with pytest.raises(MissingInputError, match="mmseqs binary"):
        cluster_sequences(input_faa, tmp_path / "clusters")
