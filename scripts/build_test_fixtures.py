"""Build deterministic test fixtures for GasRegNet."""

from __future__ import annotations

import json
import sqlite3
from pathlib import Path

FIXTURES_DIR = Path(__file__).resolve().parents[1] / "tests" / "fixtures"
MINI_EFI = FIXTURES_DIR / "mini_efi.sqlite"


def _json(values: list[str]) -> str:
    return json.dumps(values, separators=(",", ":"))


def build_mini_efi(path: Path = MINI_EFI) -> Path:
    """Build a deterministic normalized EFI-GNT SQLite fixture."""

    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        path.unlink()

    with sqlite3.connect(path) as connection:
        connection.executescript(
            """
            create table neighborhoods (
                locus_key text primary key,
                cluster_id integer not null,
                anchor_accession text not null,
                anchor_family text not null,
                organism text not null,
                taxon_id integer not null,
                contig_id text not null,
                window_size integer not null,
                is_boundary_truncated integer not null,
                marker_genes_present_json text not null,
                accessory_genes_present_json text not null
            );

            create table genes (
                locus_key text not null,
                gene_accession text not null,
                gene_order integer not null,
                start_nt integer not null,
                stop_nt integer not null,
                strand text not null,
                product_description text not null,
                pfam_ids_json text not null,
                pfam_descriptions_json text not null,
                interpro_ids_json text not null,
                interpro_descriptions_json text not null,
                functional_class text not null,
                regulator_class text not null,
                sensory_domains_json text not null,
                is_regulator_candidate integer not null,
                primary key (locus_key, gene_accession),
                foreign key (locus_key) references neighborhoods(locus_key)
            );
            """,
        )
        connection.execute(
            """
            insert into neighborhoods values
            (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                "locus_a",
                7,
                "COX_L_ANCHOR",
                "coxL",
                "Rhodospirillum rubrum",
                1085,
                "contig_1",
                5,
                0,
                _json(["coxL"]),
                _json(["coxM", "coxS"]),
            ),
        )
        connection.executemany(
            """
            insert into genes values
            (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            [
                (
                    "locus_a",
                    "REG_UP",
                    1,
                    100,
                    450,
                    "+",
                    "PAS domain transcriptional regulator",
                    _json(["PF00989", "PF00325"]),
                    _json(["PAS", "HTH"]),
                    _json(["IPR000014"]),
                    _json(["PAS-like sensor"]),
                    "regulator",
                    "one_component",
                    _json(["PAS"]),
                    1,
                ),
                (
                    "locus_a",
                    "COX_M",
                    2,
                    600,
                    980,
                    "+",
                    "carbon monoxide dehydrogenase medium subunit",
                    _json(["PF03450"]),
                    _json(["coxM"]),
                    _json([]),
                    _json([]),
                    "metabolic",
                    "none",
                    _json([]),
                    0,
                ),
                (
                    "locus_a",
                    "COX_L_ANCHOR",
                    3,
                    1100,
                    2300,
                    "+",
                    "carbon monoxide dehydrogenase large subunit",
                    _json(["PF02738"]),
                    _json(["coxL"]),
                    _json([]),
                    _json([]),
                    "anchor",
                    "none",
                    _json([]),
                    0,
                ),
                (
                    "locus_a",
                    "COX_S",
                    4,
                    2400,
                    2700,
                    "+",
                    "carbon monoxide dehydrogenase small subunit",
                    _json(["PF01799"]),
                    _json(["coxS"]),
                    _json([]),
                    _json([]),
                    "metabolic",
                    "none",
                    _json([]),
                    0,
                ),
            ],
        )
    return path


def main() -> None:
    """Build all deterministic test fixtures."""

    build_mini_efi()


if __name__ == "__main__":
    main()
