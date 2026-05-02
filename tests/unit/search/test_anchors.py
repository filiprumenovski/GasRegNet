from __future__ import annotations

from pathlib import Path

import duckdb
import polars as pl
import pytest

from gasregnet.config import load_config
from gasregnet.search.anchors import (
    _hits_for_family,
    _passes_family_guard,
    _seed_sequences,
    detect_anchors_profile,
)


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
                ('NP_TRUE.1', 'cytochrome bd oxidase subunit I', 'MLDIV', 5),
                ('NP_FALSE.1', 'hypothetical protein', 'MLDIV', 5)
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
                ('NP_TRUE.1', 'b0001', 'cydA', 'cytochrome bd oxidase subunit I'),
                ('NP_FALSE.1', 'b0002', 'foo', 'hypothetical protein')
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
    profile_dir = tmp_path / "profiles"
    profile_dir.mkdir()
    (profile_dir / "cydA.hmm").write_text("HMMER3/f\n", encoding="utf-8")

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
        profile_dir=profile_dir,
        e_value_threshold=1e-20,
        back_confirm_coverage=0.0,
    )

    assert hits.height == 1
    assert hits["dataset_name"].item() == "mini"
    assert hits["analyte"].item() == "cyd_control"
    assert hits["anchor_family"].item() == "cydA"
    assert hits["protein_accession"].item() == "NP_TRUE.1"
    assert hits["locus_tag"].item() == "b0001"
    assert hits["evidence_type"].item() == "seed_back_confirmed"
    assert hits["identity"].item() > 0.0
    assert hits["coverage"].item() > 0.0


def test_hits_for_family_filters_family_guard_vectorized(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _catalog(tmp_path)

    def fake_hmmsearch(
        _profile_hmm: Path,
        _sequences_faa: Path,
        *,
        e_value: float,
        cpus: int = 0,
    ) -> pl.DataFrame:
        return pl.DataFrame(
            {
                "query_id": ["cydA", "cydA"],
                "target_id": ["NP_TRUE.1", "NP_FALSE.1"],
                "evalue": [min(e_value, 1e-30), min(e_value, 1e-30)],
                "score": [80.0, 90.0],
                "bias": [0.0, 0.0],
                "included": [True, True],
            },
        )

    monkeypatch.setattr("gasregnet.search.anchors.hmmsearch", fake_hmmsearch)

    hits = _hits_for_family(
        dataset_name="mini",
        db=tmp_path / "mini.duckdb",
        protein_faa=tmp_path / "proteins.faa",
        analyte="CN",
        anchor_family="cydA",
        profile_hmm=tmp_path / "profiles" / "cydA.hmm",
        e_value_threshold=1e-20,
        bitscore_threshold=None,
        seeds_by_analyte={},
        back_confirm_identity=0.25,
        back_confirm_coverage=0.5,
    )

    assert hits["protein_accession"].to_list() == ["NP_TRUE.1"]
    assert hits["locus_tag"].to_list() == ["b0001"]


def test_co_seed_sequences_split_coxl_and_coos_without_downloads() -> None:
    config = load_config(Path("configs"))

    co_seeds = _seed_sequences(config)["CO"]

    assert co_seeds["coxL"].startswith("MNIQTTV")
    assert co_seeds["cooS"].startswith("MTHHDCA")
    assert co_seeds["coxL"] != co_seeds["cooS"]


def test_co_family_guard_keeps_coxl_and_coos_separate() -> None:
    coxl_row = {
        "protein_accession": "NP_COXL.1",
        "locus_tag": "coxL",
        "gene": "coxL",
        "product": "aerobic carbon monoxide dehydrogenase large subunit",
    }
    coos_row = {
        "protein_accession": "NP_COOS.1",
        "locus_tag": "cooS",
        "gene": "cooS",
        "product": "anaerobic carbon monoxide dehydrogenase catalytic subunit",
    }

    assert _passes_family_guard(coxl_row, "coxL")
    assert not _passes_family_guard(coxl_row, "cooS")
    assert _passes_family_guard(coos_row, "cooS")
    assert not _passes_family_guard(coos_row, "coxL")
