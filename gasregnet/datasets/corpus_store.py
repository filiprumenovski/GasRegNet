"""Partitioned-Parquet corpus store helpers."""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path

import duckdb
import polars as pl
import pyarrow.dataset as ds  # type: ignore[import-untyped]

TAXONOMY_COLUMNS = [
    "taxon_id",
    "superkingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
]


def _with_dataset_taxonomy(
    frame: pl.DataFrame,
    *,
    dataset_name: str,
    taxonomy: Mapping[str, object],
) -> pl.DataFrame:
    columns: list[pl.Expr] = [pl.lit(dataset_name).alias("dataset_name")]
    for column in TAXONOMY_COLUMNS:
        value = taxonomy.get(column, 0 if column == "taxon_id" else "")
        dtype = pl.Int64 if column == "taxon_id" else pl.Utf8
        columns.append(pl.lit(value).cast(dtype).alias(column))
    return frame.with_columns(columns)


def _write_dataset(
    frame: pl.DataFrame,
    root: Path,
    *,
    partition_cols: list[str],
) -> None:
    root.mkdir(parents=True, exist_ok=True)
    table = frame.to_arrow()
    ds.write_dataset(
        table,
        base_dir=str(root),
        format="parquet",
        partitioning=partition_cols,
        partitioning_flavor="hive",
        existing_data_behavior="delete_matching",
    )


def write_catalog(
    *,
    store_root: Path,
    dataset_name: str,
    proteins: pl.DataFrame,
    features: pl.DataFrame,
    taxonomy: Mapping[str, object],
    metadata: dict[str, object] | None = None,
) -> Path:
    """Write one indexed catalog into the partitioned corpus store."""

    normalized_taxonomy = {
        "taxon_id": int(str(taxonomy.get("taxon_id", 0) or 0)),
        "superkingdom": str(taxonomy.get("superkingdom", "")),
        "phylum": str(taxonomy.get("phylum", "")),
        "class": str(taxonomy.get("class", "")),
        "order": str(taxonomy.get("order", "")),
        "family": str(taxonomy.get("family", "")),
        "genus": str(taxonomy.get("genus", "")),
        "species": str(taxonomy.get("species", "")),
    }
    protein_frame = _with_dataset_taxonomy(
        proteins,
        dataset_name=dataset_name,
        taxonomy=normalized_taxonomy,
    )
    feature_frame = _with_dataset_taxonomy(
        features,
        dataset_name=dataset_name,
        taxonomy=normalized_taxonomy,
    )
    partition_cols = ["phylum", "dataset_name"]
    _write_dataset(
        protein_frame,
        store_root / "proteins",
        partition_cols=partition_cols,
    )
    _write_dataset(
        feature_frame,
        store_root / "features",
        partition_cols=partition_cols,
    )

    dataset_row: dict[str, object] = {
        "dataset_name": dataset_name,
        **normalized_taxonomy,
    }
    if metadata:
        dataset_row.update({key: str(value) for key, value in metadata.items()})
    datasets_path = store_root / "datasets.parquet"
    existing = (
        pl.read_parquet(datasets_path)
        if datasets_path.exists()
        else pl.DataFrame(
            schema={key: pl.Utf8 for key in dataset_row if key != "taxon_id"},
        )
    )
    next_row = pl.DataFrame([dataset_row])
    if not existing.is_empty() and "dataset_name" in existing.columns:
        existing = existing.filter(pl.col("dataset_name") != dataset_name)
    pl.concat([existing, next_row], how="diagonal_relaxed").write_parquet(datasets_path)
    return store_root


def register_views(
    connection: duckdb.DuckDBPyConnection,
    store_root: Path,
    *,
    view_prefix: str = "",
) -> None:
    """Register DuckDB views over the corpus-store Parquet datasets."""

    protein_glob = str(store_root / "proteins" / "**" / "*.parquet").replace("'", "''")
    feature_glob = str(store_root / "features" / "**" / "*.parquet").replace("'", "''")
    datasets_path = str(store_root / "datasets.parquet").replace("'", "''")
    connection.execute(
        f"""
        CREATE OR REPLACE VIEW {view_prefix}proteins AS
        SELECT * FROM parquet_scan('{protein_glob}', hive_partitioning=true)
        """,
    )
    connection.execute(
        f"""
        CREATE OR REPLACE VIEW {view_prefix}features AS
        SELECT * FROM parquet_scan('{feature_glob}', hive_partitioning=true)
        """,
    )
    connection.execute(
        f"""
        CREATE OR REPLACE VIEW {view_prefix}datasets AS
        SELECT * FROM parquet_scan('{datasets_path}')
        """,
    )


def catalog_counts(store_root: Path) -> pl.DataFrame:
    """Return per-dataset protein and feature counts from a corpus store."""

    with duckdb.connect() as connection:
        register_views(connection, store_root)
        return connection.execute(
            """
            SELECT
                datasets.dataset_name,
                datasets.taxon_id,
                datasets.phylum,
                count(DISTINCT proteins.protein_accession) AS n_proteins,
                count(features.feature_id) AS n_features
            FROM datasets
            LEFT JOIN proteins USING (dataset_name)
            LEFT JOIN features USING (dataset_name)
            GROUP BY datasets.dataset_name, datasets.taxon_id, datasets.phylum
            ORDER BY datasets.dataset_name
            """,
        ).pl()
