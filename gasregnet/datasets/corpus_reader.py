"""Shared DuckDB reader for partitioned corpus stores."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import duckdb
import polars as pl

FEATURE_COLUMNS = [
    "seqid",
    "source",
    "feature_type",
    "start_nt",
    "end_nt",
    "strand",
    "phase",
    "feature_id",
    "parent_id",
    "locus_tag",
    "protein_accession",
    "gene",
    "product",
]


@dataclass
class CorpusStoreHandle:
    """Context-managed DuckDB handle over a Hive-partitioned corpus store."""

    store_root: Path
    connection: duckdb.DuckDBPyConnection

    def _dataset_phylum(self, dataset_name: str) -> str:
        row = self.connection.execute(
            "select phylum from datasets where dataset_name = ?",
            [dataset_name],
        ).fetchone()
        if row is None:
            raise ValueError(f"dataset not found in corpus store: {dataset_name}")
        return str(row[0])

    def _partition_path(self, table: str, dataset_name: str) -> Path:
        phylum = self._dataset_phylum(dataset_name)
        path = (
            self.store_root
            / table
            / f"phylum={phylum}"
            / f"dataset_name={dataset_name}"
            / "part-0.parquet"
        )
        if not path.exists():
            raise FileNotFoundError(path)
        return path

    @staticmethod
    def _path_sql(path: Path) -> str:
        return str(path).replace("'", "''")

    def __enter__(self) -> CorpusStoreHandle:
        return self

    def __exit__(self, *_exc: object) -> None:
        self.close()

    def close(self) -> None:
        self.connection.close()

    def datasets(self) -> pl.DataFrame:
        return self.connection.execute(
            "select * from datasets order by dataset_name",
        ).pl()

    def dataset_taxonomy(self, dataset_name: str) -> dict[str, str | int]:
        row = self.connection.execute(
            """
            select
                taxon_id, superkingdom, phylum, class, "order",
                family, genus, species
            from datasets
            where dataset_name = ?
            """,
            [dataset_name],
        ).fetchone()
        if row is None:
            raise ValueError(f"dataset not found in corpus store: {dataset_name}")
        keys = [
            "taxon_id",
            "superkingdom",
            "phylum",
            "class",
            "order",
            "family",
            "genus",
            "species",
        ]
        values = dict(zip(keys, row, strict=True))
        values["taxon_id"] = int(str(values["taxon_id"] or 0))
        return values

    def fetch_protein_metadata(
        self,
        protein_accessions: list[str],
        *,
        dataset_name: str | None = None,
    ) -> pl.DataFrame:
        if not protein_accessions:
            return pl.DataFrame(
                schema={
                    "dataset_name": pl.Utf8,
                    "protein_accession": pl.Utf8,
                    "locus_tag": pl.Utf8,
                    "gene": pl.Utf8,
                    "product": pl.Utf8,
                    "sequence": pl.Utf8,
                },
            )
        if dataset_name:
            proteins_path = self._path_sql(
                self._partition_path("proteins", dataset_name),
            )
            features_path = self._path_sql(
                self._partition_path("features", dataset_name),
            )
            return self.connection.execute(
                f"""
                select
                    ? as dataset_name,
                    proteins.protein_accession,
                    coalesce(features.locus_tag, '') as locus_tag,
                    coalesce(features.gene, '') as gene,
                    coalesce(features.product, proteins.description) as product,
                    proteins.sequence
                from parquet_scan('{proteins_path}') as proteins
                left join parquet_scan('{features_path}') as features
                  on proteins.protein_accession = features.protein_accession
                where proteins.protein_accession in (select unnest(?::varchar[]))
                """,
                [dataset_name, protein_accessions],
            ).pl()
        params: list[object] = [protein_accessions]
        if dataset_name:
            params.append(dataset_name)
        return self.connection.execute(
            """
            select
                proteins.dataset_name,
                proteins.protein_accession,
                coalesce(features.locus_tag, '') as locus_tag,
                coalesce(features.gene, '') as gene,
                coalesce(features.product, proteins.description) as product,
                proteins.sequence
            from proteins
            left join features
              on proteins.dataset_name = features.dataset_name
             and proteins.protein_accession = features.protein_accession
            where proteins.protein_accession in (select unnest(?::varchar[]))
            """,
            params,
        ).pl()

    def fetch_proteins_for_dataset(self, dataset_name: str) -> pl.DataFrame:
        proteins_path = self._path_sql(self._partition_path("proteins", dataset_name))
        features_path = self._path_sql(self._partition_path("features", dataset_name))
        return self.connection.execute(
            f"""
            select
                ? as dataset_name,
                proteins.protein_accession,
                coalesce(features.locus_tag, '') as locus_tag,
                coalesce(features.gene, '') as gene,
                coalesce(features.product, proteins.description) as product,
                proteins.sequence
            from parquet_scan('{proteins_path}') as proteins
            left join parquet_scan('{features_path}') as features
              on proteins.protein_accession = features.protein_accession
            order by proteins.protein_accession
            """,
            [dataset_name],
        ).pl()

    def fetch_cds_features_for_seqids(
        self,
        *,
        dataset_name: str,
        seqids: list[str],
    ) -> pl.DataFrame:
        if not seqids:
            return pl.DataFrame(
                schema={
                    "dataset_name": pl.Utf8,
                    **{column: pl.Utf8 for column in FEATURE_COLUMNS},
                },
            )
        features_path = self._path_sql(self._partition_path("features", dataset_name))
        return self.connection.execute(
            f"""
            select
                ? as dataset_name,
                seqid,
                source,
                feature_type,
                start_nt,
                end_nt,
                strand,
                phase,
                coalesce(feature_id, '') as feature_id,
                coalesce(parent_id, '') as parent_id,
                coalesce(locus_tag, '') as locus_tag,
                coalesce(protein_accession, '') as protein_accession,
                coalesce(gene, '') as gene,
                coalesce(product, '') as product
            from parquet_scan('{features_path}')
            where feature_type = 'CDS'
              and seqid in (select unnest(?::varchar[]))
              and protein_accession is not null
            order by seqid, start_nt, end_nt, protein_accession
            """,
            [dataset_name, seqids],
        ).pl()

    def fetch_anchor_feature_rows(
        self,
        *,
        dataset_name: str,
        anchors: list[dict[str, object]],
    ) -> dict[int, dict[str, object]]:
        if not anchors:
            return {}
        anchor_frame = pl.DataFrame(
            [
                {
                    "anchor_row": index,
                    "protein_accession": str(anchor["protein_accession"]),
                    "locus_tag": str(anchor["locus_tag"]),
                    "gene": str(anchor["gene"]),
                }
                for index, anchor in enumerate(anchors)
            ],
            schema={
                "anchor_row": pl.Int64,
                "protein_accession": pl.Utf8,
                "locus_tag": pl.Utf8,
                "gene": pl.Utf8,
            },
        )
        self.connection.register("anchor_lookup", anchor_frame)
        features_path = self._path_sql(self._partition_path("features", dataset_name))
        frame = self.connection.execute(
            f"""
            with matched as (
                select
                    anchor_lookup.anchor_row,
                    0 as match_rank,
                    features.seqid,
                    features.source,
                    features.feature_type,
                    features.start_nt,
                    features.end_nt,
                    features.strand,
                    features.phase,
                    coalesce(features.feature_id, '') as feature_id,
                    coalesce(features.parent_id, '') as parent_id,
                    coalesce(features.locus_tag, '') as locus_tag,
                    coalesce(features.protein_accession, '') as protein_accession,
                    coalesce(features.gene, '') as gene,
                    coalesce(features.product, '') as product
                from anchor_lookup
                join parquet_scan('{features_path}') as features
                  on anchor_lookup.protein_accession <> ''
                 and features.feature_type = 'CDS'
                 and features.protein_accession = anchor_lookup.protein_accession

                union all

                select
                    anchor_lookup.anchor_row,
                    1 as match_rank,
                    features.seqid,
                    features.source,
                    features.feature_type,
                    features.start_nt,
                    features.end_nt,
                    features.strand,
                    features.phase,
                    coalesce(features.feature_id, '') as feature_id,
                    coalesce(features.parent_id, '') as parent_id,
                    coalesce(features.locus_tag, '') as locus_tag,
                    coalesce(features.protein_accession, '') as protein_accession,
                    coalesce(features.gene, '') as gene,
                    coalesce(features.product, '') as product
                from anchor_lookup
                join parquet_scan('{features_path}') as features
                  on anchor_lookup.locus_tag <> ''
                 and features.feature_type = 'CDS'
                 and features.locus_tag = anchor_lookup.locus_tag

                union all

                select
                    anchor_lookup.anchor_row,
                    2 as match_rank,
                    features.seqid,
                    features.source,
                    features.feature_type,
                    features.start_nt,
                    features.end_nt,
                    features.strand,
                    features.phase,
                    coalesce(features.feature_id, '') as feature_id,
                    coalesce(features.parent_id, '') as parent_id,
                    coalesce(features.locus_tag, '') as locus_tag,
                    coalesce(features.protein_accession, '') as protein_accession,
                    coalesce(features.gene, '') as gene,
                    coalesce(features.product, '') as product
                from anchor_lookup
                join parquet_scan('{features_path}') as features
                  on anchor_lookup.gene <> ''
                 and features.feature_type = 'CDS'
                 and features.gene = anchor_lookup.gene
            ),
            ranked as (
                select
                    *,
                    row_number() over (
                        partition by anchor_row
                        order by match_rank, start_nt, end_nt
                    ) as row_number
                from matched
            )
            select
                anchor_row,
                seqid,
                source,
                feature_type,
                start_nt,
                end_nt,
                strand,
                phase,
                feature_id,
                parent_id,
                locus_tag,
                protein_accession,
                gene,
                product
            from ranked
            where row_number = 1
            order by anchor_row
            """,
        ).pl()
        self.connection.unregister("anchor_lookup")
        return {
            int(row["anchor_row"]): {
                column: row[column] for column in FEATURE_COLUMNS
            }
            for row in frame.iter_rows(named=True)
        }

    def write_dataset_fasta(self, dataset_name: str, out_fasta: Path) -> Path:
        proteins = self.fetch_proteins_for_dataset(dataset_name)
        out_fasta.parent.mkdir(parents=True, exist_ok=True)
        with out_fasta.open("w", encoding="utf-8") as handle:
            for row in proteins.iter_rows(named=True):
                handle.write(
                    f">{row['protein_accession']} {row['product']}\n"
                    f"{row['sequence']}\n",
                )
        return out_fasta


def open_corpus_store(store_root: Path) -> CorpusStoreHandle:
    """Open a corpus store with registered DuckDB views."""

    connection = duckdb.connect()
    datasets_path = str(store_root / "datasets.parquet").replace("'", "''")
    connection.execute(
        f"""
        CREATE OR REPLACE VIEW datasets AS
        SELECT * FROM parquet_scan('{datasets_path}')
        """,
    )
    return CorpusStoreHandle(store_root=store_root, connection=connection)
