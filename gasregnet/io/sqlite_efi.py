"""EFI-GNT SQLite ingestion."""

from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import duckdb
import polars as pl

from gasregnet.errors import MissingInputError, SchemaError
from gasregnet.schemas import GenesSchema, LociSchema, validate

REQUIRED_TABLES = {"neighborhoods", "genes"}
NEIGHBORHOOD_COLUMNS = {
    "locus_key",
    "cluster_id",
    "anchor_accession",
    "anchor_family",
    "organism",
    "taxon_id",
    "contig_id",
    "window_size",
    "is_boundary_truncated",
    "marker_genes_present_json",
    "accessory_genes_present_json",
}
GENE_COLUMNS = {
    "locus_key",
    "gene_accession",
    "gene_order",
    "start_nt",
    "stop_nt",
    "strand",
    "product_description",
    "pfam_ids_json",
    "pfam_descriptions_json",
    "interpro_ids_json",
    "interpro_descriptions_json",
    "functional_class",
    "regulator_class",
    "sensory_domains_json",
    "is_regulator_candidate",
}
LOCI_SCHEMA: dict[str, Any] = {
    "locus_id": pl.Utf8,
    "analyte": pl.Utf8,
    "anchor_accession": pl.Utf8,
    "anchor_family": pl.Utf8,
    "organism": pl.Utf8,
    "taxon_id": pl.Int64,
    "cluster_id": pl.Int32,
    "contig_id": pl.Utf8,
    "window_size": pl.Int32,
    "is_boundary_truncated": pl.Boolean,
    "marker_genes_present": pl.List(pl.Utf8),
    "accessory_genes_present": pl.List(pl.Utf8),
    "locus_score": pl.Float64,
    "locus_confidence": pl.Utf8,
    "taxonomic_context_score": pl.Float64,
    "operon_integrity_score": pl.Float64,
    "created_at": pl.Datetime("us"),
}
GENES_SCHEMA: dict[str, Any] = {
    "locus_id": pl.Utf8,
    "gene_accession": pl.Utf8,
    "relative_index": pl.Int32,
    "relative_start": pl.Int64,
    "relative_stop": pl.Int64,
    "strand": pl.Utf8,
    "product_description": pl.Utf8,
    "pfam_ids": pl.List(pl.Utf8),
    "pfam_descriptions": pl.List(pl.Utf8),
    "interpro_ids": pl.List(pl.Utf8),
    "interpro_descriptions": pl.List(pl.Utf8),
    "functional_class": pl.Utf8,
    "regulator_class": pl.Utf8,
    "sensory_domains": pl.List(pl.Utf8),
    "is_anchor": pl.Boolean,
    "is_regulator_candidate": pl.Boolean,
}


def _connect(path: Path) -> duckdb.DuckDBPyConnection:
    if not path.exists():
        raise MissingInputError(f"SQLite file does not exist: {path}")
    connection = duckdb.connect()
    quoted_path = str(path).replace("'", "''")
    try:
        connection.execute("INSTALL sqlite")
        connection.execute("LOAD sqlite")
        connection.execute(f"ATTACH '{quoted_path}' AS efi (TYPE SQLITE)")
    except duckdb.Error as exc:
        connection.close()
        raise SchemaError(
            "SQLite file is missing required table(s): genes, neighborhoods",
        ) from exc
    return connection


def _table_names(connection: duckdb.DuckDBPyConnection) -> set[str]:
    rows = connection.execute(
        """
        select distinct table_name
        from information_schema.columns
        where table_catalog = 'efi'
        """,
    ).fetchall()
    return {str(row[0]) for row in rows}


def _columns(connection: duckdb.DuckDBPyConnection, table: str) -> set[str]:
    rows = connection.execute(
        """
        select column_name
        from information_schema.columns
        where table_catalog = 'efi' and table_name = ?
        """,
        [table],
    ).fetchall()
    return {str(row[0]) for row in rows}


def _check_sqlite_shape(connection: duckdb.DuckDBPyConnection) -> None:
    tables = _table_names(connection)
    missing_tables = REQUIRED_TABLES - tables
    if missing_tables:
        missing = ", ".join(sorted(missing_tables))
        raise SchemaError(f"SQLite file is missing required table(s): {missing}")

    missing_neighborhood_columns = NEIGHBORHOOD_COLUMNS - _columns(
        connection,
        "neighborhoods",
    )
    if missing_neighborhood_columns:
        missing = ", ".join(sorted(missing_neighborhood_columns))
        raise SchemaError(f"neighborhoods table is missing column(s): {missing}")

    missing_gene_columns = GENE_COLUMNS - _columns(connection, "genes")
    if missing_gene_columns:
        missing = ", ".join(sorted(missing_gene_columns))
        raise SchemaError(f"genes table is missing column(s): {missing}")


def _json_list(value: object) -> list[str]:
    if value is None:
        return []
    if isinstance(value, str):
        parsed = json.loads(value)
        if not isinstance(parsed, list) or not all(
            isinstance(item, str) for item in parsed
        ):
            raise SchemaError(f"expected JSON list of strings, got {value!r}")
        return parsed
    raise SchemaError(f"expected JSON list string, got {type(value).__name__}")


def _fetch_dicts(
    connection: duckdb.DuckDBPyConnection,
    query: str,
    parameters: tuple[Any, ...] = (),
) -> list[dict[str, Any]]:
    result = connection.execute(query, parameters)
    if result.description is None:
        return []
    columns = [column[0] for column in result.description]
    rows = result.fetchall()
    return [dict(zip(columns, row, strict=True)) for row in rows]


def _cluster_clause(cluster_filter: list[int] | None) -> tuple[str, tuple[Any, ...]]:
    if cluster_filter is None:
        return "", ()
    if not cluster_filter:
        return " where 0", ()
    placeholders = ", ".join("?" for _ in cluster_filter)
    return f" where cluster_id in ({placeholders})", tuple(cluster_filter)


def _build_loci(neighborhoods: list[dict[str, Any]], analyte: str) -> pl.DataFrame:
    if not neighborhoods:
        return pl.DataFrame(schema=LOCI_SCHEMA)

    now = datetime.now(UTC).replace(tzinfo=None)
    rows = [
        {
            "locus_id": f"{analyte}_{row['cluster_id']}_{row['anchor_accession']}",
            "analyte": analyte,
            "anchor_accession": row["anchor_accession"],
            "anchor_family": row["anchor_family"],
            "organism": row["organism"],
            "taxon_id": row["taxon_id"],
            "cluster_id": row["cluster_id"],
            "contig_id": row["contig_id"],
            "window_size": row["window_size"],
            "is_boundary_truncated": bool(row["is_boundary_truncated"]),
            "marker_genes_present": _json_list(row["marker_genes_present_json"]),
            "accessory_genes_present": _json_list(
                row["accessory_genes_present_json"],
            ),
            "locus_score": 0.0,
            "locus_confidence": "low",
            "taxonomic_context_score": 0.0,
            "operon_integrity_score": 0.0,
            "created_at": now,
        }
        for row in neighborhoods
    ]
    return pl.DataFrame(
        rows,
        schema_overrides={
            "cluster_id": pl.Int32,
            "window_size": pl.Int32,
            "created_at": pl.Datetime("us"),
        },
    )


def _anchor_lookup(neighborhoods: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    return {str(row["locus_key"]): row for row in neighborhoods}


def _build_genes(
    neighborhoods: list[dict[str, Any]],
    genes: list[dict[str, Any]],
    analyte: str,
) -> pl.DataFrame:
    if not neighborhoods:
        return pl.DataFrame(schema=GENES_SCHEMA)

    anchors = _anchor_lookup(neighborhoods)
    genes_by_locus: dict[str, list[dict[str, Any]]] = {}
    for gene in genes:
        genes_by_locus.setdefault(str(gene["locus_key"]), []).append(gene)

    rows: list[dict[str, Any]] = []
    for locus_key, locus_genes in genes_by_locus.items():
        neighborhood = anchors[locus_key]
        anchor_accession = str(neighborhood["anchor_accession"])
        sorted_genes = sorted(locus_genes, key=lambda row: int(row["gene_order"]))
        anchor_rows = [
            row for row in sorted_genes if row["gene_accession"] == anchor_accession
        ]
        if len(anchor_rows) != 1:
            raise SchemaError(
                f"locus {locus_key} must contain exactly one anchor gene "
                f"{anchor_accession!r}",
            )
        anchor_row = anchor_rows[0]
        anchor_order = int(anchor_row["gene_order"])
        anchor_start = int(anchor_row["start_nt"])
        locus_id = f"{analyte}_{neighborhood['cluster_id']}_{anchor_accession}"

        for gene in sorted_genes:
            is_anchor = gene["gene_accession"] == anchor_accession
            rows.append(
                {
                    "locus_id": locus_id,
                    "gene_accession": gene["gene_accession"],
                    "relative_index": int(gene["gene_order"]) - anchor_order,
                    "relative_start": int(gene["start_nt"]) - anchor_start,
                    "relative_stop": int(gene["stop_nt"]) - anchor_start,
                    "strand": gene["strand"],
                    "product_description": gene["product_description"],
                    "pfam_ids": _json_list(gene["pfam_ids_json"]),
                    "pfam_descriptions": _json_list(gene["pfam_descriptions_json"]),
                    "interpro_ids": _json_list(gene["interpro_ids_json"]),
                    "interpro_descriptions": _json_list(
                        gene["interpro_descriptions_json"],
                    ),
                    "functional_class": gene["functional_class"],
                    "regulator_class": gene["regulator_class"],
                    "sensory_domains": _json_list(gene["sensory_domains_json"]),
                    "is_anchor": is_anchor,
                    "is_regulator_candidate": bool(gene["is_regulator_candidate"]),
                },
            )

    return pl.DataFrame(
        rows,
        schema_overrides={
            "relative_index": pl.Int32,
            "sensory_domains": pl.List(pl.Utf8),
        },
    )


def read_efi_sqlite(
    path: Path,
    analyte: str,
    cluster_filter: list[int] | None = None,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Read an EFI-GNT SQLite export into validated loci and genes tables."""

    with _connect(path) as connection:
        _check_sqlite_shape(connection)
        clause, parameters = _cluster_clause(cluster_filter)
        neighborhoods = _fetch_dicts(
            connection,
            f"select * from efi.neighborhoods{clause} order by cluster_id, locus_key",
            parameters,
        )
        genes = _fetch_dicts(
            connection,
            """
            select genes.*
            from efi.genes
            join efi.neighborhoods using (locus_key)
            """
            f"{clause}"
            " order by genes.locus_key, genes.gene_order",
            parameters,
        )

    loci_df = _build_loci(neighborhoods, analyte)
    genes_df = _build_genes(neighborhoods, genes, analyte)
    return validate(loci_df, LociSchema), validate(genes_df, GenesSchema)
