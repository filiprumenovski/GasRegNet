from __future__ import annotations

from pathlib import Path

import duckdb

from gasregnet.annotation.ncbi_taxonomy import build_taxonomy_db, lookup_lineage

FIXTURE_TAXDUMP = Path("tests/fixtures/taxdump_minimal")


def test_build_taxonomy_db_writes_lineage_tables(tmp_path: Path) -> None:
    db = build_taxonomy_db(FIXTURE_TAXDUMP, tmp_path / "taxonomy.duckdb")

    assert db == tmp_path / "taxonomy.duckdb"
    assert db.exists()

    with duckdb.connect(str(db), read_only=True) as connection:
        tables = {
            row[0]
            for row in connection.execute("SHOW TABLES").fetchall()
        }
        lineage_rows = connection.execute(
            """
            SELECT ancestor_rank, ancestor_name, distance
            FROM lineage
            WHERE taxon_id = 1085
            ORDER BY distance
            """,
        ).fetchall()
        pivot = connection.execute(
            """
            SELECT scientific_name, phylum, class, genus, species
            FROM lineage_pivot
            WHERE taxon_id = 1085
            """,
        ).fetchone()

    assert tables == {"lineage", "lineage_pivot"}
    assert lineage_rows[:3] == [
        ("species", "Rhodospirillum rubrum", 0),
        ("genus", "Rhodospirillum", 1),
        ("family", "Rhodospirillaceae", 2),
    ]
    assert pivot == (
        "Rhodospirillum rubrum",
        "Pseudomonadota",
        "Alphaproteobacteria",
        "Rhodospirillum",
        "Rhodospirillum rubrum",
    )


def test_lookup_lineage_returns_populated_rank_names(tmp_path: Path) -> None:
    db = build_taxonomy_db(FIXTURE_TAXDUMP, tmp_path / "taxonomy.duckdb")

    lineage = lookup_lineage(db, 1085)

    assert lineage == {
        "scientific_name": "Rhodospirillum rubrum",
        "superkingdom": "Bacteria",
        "phylum": "Pseudomonadota",
        "class": "Alphaproteobacteria",
        "order": "Rhodospirillales",
        "family": "Rhodospirillaceae",
        "genus": "Rhodospirillum",
        "species": "Rhodospirillum rubrum",
    }


def test_lookup_lineage_returns_empty_dict_for_unknown_taxon(tmp_path: Path) -> None:
    db = build_taxonomy_db(FIXTURE_TAXDUMP, tmp_path / "taxonomy.duckdb")

    assert lookup_lineage(db, 999999) == {}
