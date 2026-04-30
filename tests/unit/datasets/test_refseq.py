from __future__ import annotations

import gzip
from pathlib import Path

import duckdb
import polars as pl
import pytest

from gasregnet.datasets.refseq import (
    detect_refseq_anchor_hits,
    extract_refseq_neighborhoods,
    index_refseq_corpus,
    index_refseq_dataset,
    query_refseq_catalog,
    query_refseq_corpus,
    scan_refseq_corpus,
    summarize_refseq_corpus,
)


def test_index_refseq_dataset_writes_duckdb_catalog(tmp_path: Path) -> None:
    protein_faa = tmp_path / "proteins.faa.gz"
    with gzip.open(protein_faa, "wt", encoding="utf-8") as handle:
        handle.write(
            ">NP_000001.1 carbon monoxide dehydrogenase [Example bacterium]\n" "MAGA\n",
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


def test_refseq_corpus_manifest_indexes_summarizes_and_queries(tmp_path: Path) -> None:
    protein_faa = tmp_path / "proteins.faa.gz"
    with gzip.open(protein_faa, "wt", encoding="utf-8") as handle:
        handle.write(">NP_000002.1 cytochrome bd [Example bacterium]\nMAGATA\n")
    gff = tmp_path / "genome.gff.gz"
    with gzip.open(gff, "wt", encoding="utf-8") as handle:
        handle.write(
            "##gff-version 3\n"
            "ctg\tRefSeq\tCDS\t10\t99\t.\t+\t0\t"
            "ID=cds-NP_000002.1;locus_tag=GRN_0002;protein_id=NP_000002.1;"
            "gene=cydA;product=cytochrome%20bd-I%20subunit%201\n",
        )
    manifest = tmp_path / "catalogs.yaml"
    manifest.write_text(
        """
catalogs:
  - dataset_name: mini
    protein_faa: proteins.faa.gz
    gff: genome.gff.gz
    out_db: mini.duckdb
""",
        encoding="utf-8",
    )

    outputs = index_refseq_corpus(manifest, root=tmp_path)
    summary = summarize_refseq_corpus(manifest, root=tmp_path)
    results = query_refseq_corpus(manifest, "cydA", root=tmp_path)
    scan_config = tmp_path / "scan.yaml"
    scan_config.write_text(
        """
targets:
  - analyte: CN
    terms:
      - cydA
""",
        encoding="utf-8",
    )
    scan = scan_refseq_corpus(manifest, scan_config, root=tmp_path)

    assert outputs == [tmp_path / "mini.duckdb"]
    assert summary["n_proteins"].item() == 1
    assert summary["n_features_with_protein"].item() == 1
    assert results["dataset_name"].item() == "mini"
    assert results["protein_accession"].item() == "NP_000002.1"
    assert scan["analyte"].item() == "CN"
    assert scan["gene"].item() == "cydA"


def test_detect_anchors_and_extract_refseq_neighborhoods(tmp_path: Path) -> None:
    protein_faa = tmp_path / "proteins.faa.gz"
    with gzip.open(protein_faa, "wt", encoding="utf-8") as handle:
        handle.write(
            ">NP_000010.1 upstream protein\nMA\n"
            ">NP_000011.1 cytochrome bd subunit\nMAGATA\n"
            ">NP_000012.1 downstream protein\nMAA\n",
        )
    gff = tmp_path / "genome.gff.gz"
    with gzip.open(gff, "wt", encoding="utf-8") as handle:
        handle.write(
            "##gff-version 3\n"
            "ctg\tRefSeq\tCDS\t10\t99\t.\t+\t0\t"
            "ID=cds-NP_000010.1;locus_tag=GRN_0010;protein_id=NP_000010.1;"
            "gene=up;product=upstream%20protein\n"
            "ctg\tRefSeq\tCDS\t110\t199\t.\t+\t0\t"
            "ID=cds-NP_000011.1;locus_tag=GRN_0011;protein_id=NP_000011.1;"
            "gene=cydA;product=cytochrome%20bd-I%20subunit%201\n"
            "ctg\tRefSeq\tCDS\t210\t299\t.\t-\t0\t"
            "ID=cds-NP_000012.1;locus_tag=GRN_0012;protein_id=NP_000012.1;"
            "gene=down;product=downstream%20protein\n",
        )
    manifest = tmp_path / "catalogs.yaml"
    manifest.write_text(
        """
catalogs:
  - dataset_name: mini
    protein_faa: proteins.faa.gz
    gff: genome.gff.gz
    out_db: mini.duckdb
""",
        encoding="utf-8",
    )
    scan_config = tmp_path / "scan.yaml"
    scan_config.write_text(
        """
targets:
  - analyte: CN
    terms:
      - cydA
""",
        encoding="utf-8",
    )
    index_refseq_corpus(manifest, root=tmp_path)

    anchor_hits = detect_refseq_anchor_hits(
        manifest,
        scan_config,
        root=tmp_path,
        mode="smoke",
    )
    loci, genes = extract_refseq_neighborhoods(
        anchor_hits,
        manifest,
        root=tmp_path,
        window_genes=1,
    )

    assert anchor_hits["anchor_family"].item() == "cydA"
    assert anchor_hits["evidence_type"].item() == "term_scan"
    assert loci.height == 1
    assert loci["anchor_accession"].item() == "NP_000011.1"
    assert genes.height == 3
    assert genes.filter(genes["is_anchor"])["gene_accession"].item() == "NP_000011.1"
    assert genes["relative_index"].to_list() == [-1, 0, 1]


def test_extract_refseq_neighborhoods_batches_catalog_connections(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    protein_faa = tmp_path / "proteins.faa.gz"
    with gzip.open(protein_faa, "wt", encoding="utf-8") as handle:
        handle.write(
            ">NP_000020.1 cytochrome bd subunit\nMAGATA\n"
            ">NP_000021.1 carbon monoxide dehydrogenase\nMAGATA\n",
        )
    gff = tmp_path / "genome.gff.gz"
    with gzip.open(gff, "wt", encoding="utf-8") as handle:
        handle.write(
            "##gff-version 3\n"
            "ctg\tRefSeq\tCDS\t10\t99\t.\t+\t0\t"
            "ID=cds-NP_000020.1;locus_tag=GRN_0020;protein_id=NP_000020.1;"
            "gene=cydA;product=cytochrome%20bd-I%20subunit%201\n"
            "ctg\tRefSeq\tCDS\t110\t199\t.\t+\t0\t"
            "ID=cds-NP_000021.1;locus_tag=GRN_0021;protein_id=NP_000021.1;"
            "gene=coxL;product=carbon%20monoxide%20dehydrogenase\n",
        )
    manifest = tmp_path / "catalogs.yaml"
    manifest.write_text(
        """
catalogs:
  - dataset_name: mini
    protein_faa: proteins.faa.gz
    gff: genome.gff.gz
    out_db: mini.duckdb
""",
        encoding="utf-8",
    )
    index_refseq_corpus(manifest, root=tmp_path)
    anchor_hits = pl.DataFrame(
        data={
            "dataset_name": ["mini", "mini"],
            "analyte": ["CN", "CO"],
            "anchor_family": ["cydA", "coxL"],
            "protein_accession": ["NP_000020.1", "NP_000021.1"],
            "locus_tag": ["GRN_0020", "GRN_0021"],
            "gene": ["cydA", "coxL"],
            "product": [
                "cytochrome bd-I subunit 1",
                "carbon monoxide dehydrogenase",
            ],
            "bitscore": [None, None],
            "e_value": [None, None],
            "identity": [None, None],
            "coverage": [None, None],
            "evidence_type": ["term_scan", "term_scan"],
        },
        schema_overrides={
            "bitscore": pl.Float64,
            "e_value": pl.Float64,
            "identity": pl.Float64,
            "coverage": pl.Float64,
        },
    )
    original_connect = duckdb.connect
    connect_calls = 0

    def counting_connect(*args: object, **kwargs: object) -> duckdb.DuckDBPyConnection:
        nonlocal connect_calls
        connect_calls += 1
        return original_connect(*args, **kwargs)

    monkeypatch.setattr("gasregnet.datasets.refseq.duckdb.connect", counting_connect)

    loci, genes = extract_refseq_neighborhoods(
        anchor_hits,
        manifest,
        root=tmp_path,
        window_genes=0,
    )

    assert connect_calls == 1
    assert loci.height == 2
    assert genes.height == 2
    assert genes["relative_index"].to_list() == [0, 0]
