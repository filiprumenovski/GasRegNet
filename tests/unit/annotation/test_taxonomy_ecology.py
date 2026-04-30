from __future__ import annotations

import polars as pl

from gasregnet.annotation.ecology import score_taxonomic_context
from gasregnet.annotation.taxonomy import annotate_taxonomy
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
