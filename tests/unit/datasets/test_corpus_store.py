from __future__ import annotations

import duckdb
import polars as pl

from gasregnet.datasets.corpus_store import (
    catalog_counts,
    register_views,
    write_catalog,
)


def test_write_catalog_creates_partitioned_parquet_store(tmp_path):
    store = tmp_path / "store"
    proteins = pl.DataFrame(
        {
            "protein_accession": ["NP_1"],
            "description": ["carbon monoxide dehydrogenase"],
            "sequence": ["MAGA"],
            "length_aa": [4],
        },
    )
    features = pl.DataFrame(
        {
            "seqid": ["ctg"],
            "source": ["RefSeq"],
            "feature_type": ["CDS"],
            "start_nt": [1],
            "end_nt": [12],
            "strand": ["+"],
            "phase": [0],
            "feature_id": ["cds-NP_1"],
            "parent_id": ["gene-coxL"],
            "locus_tag": ["GRN_1"],
            "protein_accession": ["NP_1"],
            "gene": ["coxL"],
            "product": ["CO dehydrogenase"],
        },
    )

    write_catalog(
        store_root=store,
        dataset_name="mini",
        proteins=proteins,
        features=features,
        taxonomy={
            "taxon_id": 511145,
            "superkingdom": "Bacteria",
            "phylum": "Pseudomonadota",
            "class": "Gammaproteobacteria",
            "order": "Enterobacterales",
            "family": "Enterobacteriaceae",
            "genus": "Escherichia",
            "species": "Escherichia coli",
        },
    )

    assert list((store / "proteins").glob("phylum=*/dataset_name=*/*.parquet"))
    assert list((store / "features").glob("phylum=*/dataset_name=*/*.parquet"))
    counts = catalog_counts(store)
    assert counts.select("dataset_name", "n_proteins", "n_features").row(0) == (
        "mini",
        1,
        1,
    )
    with duckdb.connect() as connection:
        register_views(connection, store)
        row = connection.execute(
            "select taxon_id, phylum, sequence from proteins",
        ).fetchone()
    assert row == (511145, "Pseudomonadota", "MAGA")
