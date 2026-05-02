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


def test_score_taxonomic_context_matches_mixed_organism_and_taxon_rows() -> None:
    loci = pl.concat(
        [
            loci_frame().with_columns(
                pl.lit("known-by-organism").alias("locus_id"),
                pl.lit("Matched organism").alias("organism"),
                pl.lit(111).cast(pl.Int64).alias("taxon_id"),
            ),
            loci_frame().with_columns(
                pl.lit("known-by-taxon").alias("locus_id"),
                pl.lit("Different organism").alias("organism"),
                pl.lit(222).cast(pl.Int64).alias("taxon_id"),
            ),
            loci_frame().with_columns(
                pl.lit("unknown").alias("locus_id"),
                pl.lit("Unknown organism").alias("organism"),
                pl.lit(333).cast(pl.Int64).alias("taxon_id"),
            ),
        ],
        how="vertical",
    )
    known = pl.DataFrame(
        {
            "organism": ["Matched organism", "Other organism"],
            "taxon_id": [999, 222],
            "evidence": ["literature", "literature"],
            "pmid": ["00000000", "00000000"],
            "notes": ["", ""],
        },
    )

    scored = score_taxonomic_context(loci, known, matched_score=0.5).sort("locus_id")

    assert scored["taxonomic_context_score"].to_list() == [0.5, 0.5, 0.0]


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


def test_score_taxonomic_context_by_analyte_scores_with_vectorized_joins(
    tmp_path: Path,
) -> None:
    co_known = tmp_path / "co_known.csv"
    no_known = tmp_path / "no_known.csv"
    co_known.write_text(
        "organism,taxon_id,evidence,pmid,notes\n"
        "CO organism,111,literature,00000000,\n",
        encoding="utf-8",
    )
    no_known.write_text(
        "organism,taxon_id,evidence,pmid,notes\n"
        "Other organism,222,literature,00000000,\n",
        encoding="utf-8",
    )
    config = load_config(Path("configs"))
    analytes = [
        config.analytes[0].model_copy(
            update={"analyte": "CO", "known_organisms_table": co_known},
        ),
        config.analytes[1].model_copy(
            update={"analyte": "NO", "known_organisms_table": no_known},
        ),
    ]
    loci = pl.concat(
        [
            loci_frame().with_columns(
                pl.lit("no-match-by-taxon").alias("locus_id"),
                pl.lit("NO").alias("analyte"),
                pl.lit("Different organism").alias("organism"),
                pl.lit(222).cast(pl.Int64).alias("taxon_id"),
            ),
            loci_frame().with_columns(
                pl.lit("co-match-by-organism").alias("locus_id"),
                pl.lit("CO").alias("analyte"),
                pl.lit("CO organism").alias("organism"),
                pl.lit(999).cast(pl.Int64).alias("taxon_id"),
            ),
            loci_frame().with_columns(
                pl.lit("o2-unconfigured").alias("locus_id"),
                pl.lit("O2").alias("analyte"),
                pl.lit("CO organism").alias("organism"),
                pl.lit(111).cast(pl.Int64).alias("taxon_id"),
            ),
        ],
        how="vertical",
    )

    scored = score_taxonomic_context_by_analyte(
        loci,
        analytes,
        matched_score=0.25,
    )

    assert scored["locus_id"].to_list() == [
        "co-match-by-organism",
        "no-match-by-taxon",
        "o2-unconfigured",
    ]
    assert scored["taxonomic_context_score"].to_list() == [0.25, 0.25, 0.0]
