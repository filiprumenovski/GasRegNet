"""Taxonomy annotation helpers."""

from __future__ import annotations

import polars as pl

TAXONOMY_COLUMNS = ["phylum", "class", "order", "family", "genus"]


def annotate_taxonomy(loci: pl.DataFrame, taxonomy_table: pl.DataFrame) -> pl.DataFrame:
    """Join taxonomic lineage columns onto loci by NCBI taxon ID."""

    required = {"taxon_id", *TAXONOMY_COLUMNS}
    missing = required - set(taxonomy_table.columns)
    if missing:
        missing_text = ", ".join(sorted(missing))
        msg = f"taxonomy table is missing column(s): {missing_text}"
        raise ValueError(msg)

    lineage = taxonomy_table.select(["taxon_id", *TAXONOMY_COLUMNS]).unique(
        subset=["taxon_id"],
        keep="first",
    )
    return loci.join(lineage, on="taxon_id", how="left")
