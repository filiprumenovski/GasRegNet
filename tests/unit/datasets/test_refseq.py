from __future__ import annotations

import gzip
from pathlib import Path

import duckdb

from gasregnet.datasets.refseq import index_refseq_dataset, query_refseq_catalog


def test_index_refseq_dataset_writes_duckdb_catalog(tmp_path: Path) -> None:
    protein_faa = tmp_path / "proteins.faa.gz"
    with gzip.open(protein_faa, "wt", encoding="utf-8") as handle:
        handle.write(
            ">NP_000001.1 carbon monoxide dehydrogenase [Example bacterium]\n"
            "MAGA\n",
        )
    gff = tmp_path / "genome.gff.gz"
    with gzip.open(gff, "wt", encoding="utf-8") as handle:
        handle.write(
            "##gff-version 3\n"
            "ctg\tRefSeq\tCDS\t10\t99\t.\t+\t0\t"
            "ID=cds-NP_000001.1;Parent=gene-coxL;locus_tag=GRN_0001;"
            "protein_id=NP_000001.1;gene=coxL;product=CO%20dehydrogenase\n",
        )
    out_db = tmp_path / "refseq.duckdb"

    index_refseq_dataset(
        protein_faa=protein_faa,
        gff=gff,
        out_db=out_db,
        dataset_name="mini_refseq",
    )

    with duckdb.connect(str(out_db)) as connection:
        assert connection.execute("select count(*) from proteins").fetchone() == (1,)
        assert connection.execute("select count(*) from features").fetchone() == (1,)
        row = connection.execute(
            """
            select proteins.protein_accession, features.locus_tag, features.product
            from proteins
            join features using (protein_accession)
            """,
        ).fetchone()
        assert row == ("NP_000001.1", "GRN_0001", "CO dehydrogenase")
        metadata_rows = connection.execute("select key, value from metadata").fetchall()
        metadata = dict(metadata_rows)
        assert metadata["dataset_name"] == "mini_refseq"
        assert metadata["n_proteins"] == "1"

    results = query_refseq_catalog(out_db, "coxL")

    assert results.height == 1
    assert results["protein_accession"].item() == "NP_000001.1"
    assert results["gene"].item() == "coxL"
