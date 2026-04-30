from __future__ import annotations

from pathlib import Path

import polars as pl

from gasregnet.annotation.ecology import (
    score_taxonomic_context,
    score_taxonomic_context_by_analyte,
)
from gasregnet.annotation.taxonomy import annotate_taxonomy
from gasregnet.config import load_config
from tests.unit.test_schemas import loci_frame


def test_annotate_taxonomy_joins_lineage_columns() -> None:
    taxonomy = pl.DataFrame(
        {
            "taxon_id": [1085],
            "phylum": ["Pseudomonadota"],
            "class": ["Alphaproteobacteria"],
            "order": ["Rhodospirillales"],
            "family": ["Rhodospirillaceae"],
            "genus": ["Rhodospirillum"],
        },
    )

    annotated = annotate_taxonomy(loci_frame(), taxonomy)

    assert annotated["phylum"].item() == "Pseudomonadota"
    assert annotated["genus"].item() == "Rhodospirillum"


def test_score_taxonomic_context_matches_known_taxon() -> None:
    known = pl.DataFrame(
        {
            "organism": ["Different name"],
            "taxon_id": [1085],
            "evidence": ["literature"],
            "pmid": ["00000000"],
            "notes": [""],
        },
    )

    scored = score_taxonomic_context(loci_frame(), known, matched_score=0.75)

    assert scored["taxonomic_context_score"].item() == 0.75


def test_score_taxonomic_context_scores_unknown_as_zero() -> None:
    known = pl.DataFrame(
        {
            "organism": ["Different name"],
            "taxon_id": [999],
            "evidence": ["literature"],
            "pmid": ["00000000"],
            "notes": [""],
        },
    )

    scored = score_taxonomic_context(loci_frame(), known)

    assert scored["taxonomic_context_score"].item() == 0.0


def test_score_taxonomic_context_by_analyte_uses_configured_tables(
    tmp_path: Path,
) -> None:
    known = tmp_path / "known.csv"
    known.write_text(
        "organism,taxon_id,evidence,pmid,notes\n"
        "Rhodospirillum rubrum,1085,literature,00000000,\n",
        encoding="utf-8",
    )
    config = load_config(Path("configs"))
    analyte = config.analytes[0].model_copy(
        update={"known_organisms_table": known},
    )

    scored = score_taxonomic_context_by_analyte(loci_frame(), [analyte])

    assert scored["taxonomic_context_score"].item() == 1.0
