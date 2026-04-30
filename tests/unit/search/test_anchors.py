from __future__ import annotations

from pathlib import Path

import duckdb
import polars as pl
import pytest

from gasregnet.search.anchors import detect_anchors_profile


def _catalog(tmp_path: Path) -> Path:
    protein_faa = tmp_path / "proteins.faa"
    protein_faa.write_text(">NP_TRUE.1 true anchor\nMLDIV\n", encoding="utf-8")
    db = tmp_path / "mini.duckdb"
    with duckdb.connect(str(db)) as connection:
        connection.execute(
            """
            create table proteins (
                protein_accession varchar,
                description varchar,
                sequence varchar,
                length_aa bigint
            )
            """,
        )
        connection.execute(
            """
            insert into proteins values
                ('NP_TRUE.1', 'cytochrome bd oxidase subunit I', 'MLDIV', 5)
            """,
        )
        connection.execute(
            """
            create table features (
                protein_accession varchar,
                locus_tag varchar,
                gene varchar,
                product varchar
            )
            """,
        )
        connection.execute(
            """
            insert into features values
                ('NP_TRUE.1', 'b0001', 'cydA', 'cytochrome bd oxidase subunit I')
            """,
        )
    manifest = tmp_path / "catalogs.yaml"
    manifest.write_text(
        """
catalogs:
  - dataset_name: mini
    protein_faa: proteins.faa
    gff: genome.gff
    out_db: mini.duckdb
""",
        encoding="utf-8",
    )
    return manifest


def test_detect_anchors_profile_normalizes_hmmer_hits(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    manifest = _catalog(tmp_path)

    def fake_hmmsearch(
        profile_hmm: Path,
        _sequences_faa: Path,
        *,
        e_value: float,
        cpus: int = 0,
    ) -> pl.DataFrame:
        if profile_hmm.name != "cydA.hmm":
            return pl.DataFrame(
                schema={
                    "query_id": pl.Utf8,
                    "target_id": pl.Utf8,
                    "evalue": pl.Float64,
                    "score": pl.Float64,
                    "bias": pl.Float64,
                    "included": pl.Boolean,
                },
            )
        return pl.DataFrame(
            {
                "query_id": ["cydA"],
                "target_id": ["NP_TRUE.1"],
                "evalue": [min(e_value, 1e-30)],
                "score": [80.0],
                "bias": [0.0],
                "included": [True],
            },
        )

    monkeypatch.setattr("gasregnet.search.anchors.hmmsearch", fake_hmmsearch)

    hits = detect_anchors_profile(
        manifest,
        config=Path("configs"),
        root=tmp_path,
        profile_dir=tmp_path / "profiles",
        e_value_threshold=1e-20,
        back_confirm_coverage=0.0,
    )

    assert hits.height == 1
    assert hits["dataset_name"].item() == "mini"
    assert hits["analyte"].item() == "CN"
    assert hits["anchor_family"].item() == "cydA"
    assert hits["protein_accession"].item() == "NP_TRUE.1"
    assert hits["locus_tag"].item() == "b0001"
    assert hits["evidence_type"].item() == "seed_back_confirmed"
    assert hits["identity"].item() > 0.0
    assert hits["coverage"].item() > 0.0
