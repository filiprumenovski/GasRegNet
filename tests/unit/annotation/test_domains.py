from __future__ import annotations

import polars as pl

from gasregnet.annotation.domains import annotate_domains
from tests.unit.test_schemas import genes_frame


def test_annotate_domains_populates_pfam_and_interpro_lists() -> None:
    genes = genes_frame().with_columns(
        pl.lit([]).cast(pl.List(pl.Utf8)).alias("pfam_ids"),
        pl.lit([]).cast(pl.List(pl.Utf8)).alias("pfam_descriptions"),
        pl.lit([]).cast(pl.List(pl.Utf8)).alias("interpro_ids"),
        pl.lit([]).cast(pl.List(pl.Utf8)).alias("interpro_descriptions"),
    )
    pfam = pl.DataFrame(
        {
            "gene_accession": ["anchor", "anchor"],
            "pfam_id": ["PF02738", "PF03450"],
            "pfam_description": ["coxL", "coxM-like"],
        },
    )
    interpro = pl.DataFrame(
        {
            "gene_accession": ["anchor"],
            "interpro_id": ["IPR000001"],
            "interpro_description": ["molybdopterin oxidoreductase"],
        },
    )

    annotated = annotate_domains(genes, pfam, interpro)

    assert annotated["pfam_ids"].to_list()[0] == ["PF02738", "PF03450"]
    assert annotated["interpro_ids"].to_list()[0] == ["IPR000001"]


def test_annotate_domains_handles_empty_annotation_tables() -> None:
    empty_pfam = pl.DataFrame(
        schema={
            "gene_accession": pl.Utf8,
            "pfam_id": pl.Utf8,
            "pfam_description": pl.Utf8,
        },
    )
    empty_interpro = pl.DataFrame(
        schema={
            "gene_accession": pl.Utf8,
            "interpro_id": pl.Utf8,
            "interpro_description": pl.Utf8,
        },
    )

    annotated = annotate_domains(genes_frame(), empty_pfam, empty_interpro)

    assert annotated["pfam_ids"].to_list()[0] == []
    assert annotated["interpro_ids"].to_list()[0] == []
