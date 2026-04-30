"""Domain annotation helpers."""

from __future__ import annotations

from typing import Any

import polars as pl
import structlog

from gasregnet.schemas import GenesSchema, validate

log = structlog.get_logger(__name__)

GENES_SCHEMA_OVERRIDES: dict[str, Any] = {
    "relative_index": pl.Int32,
    "relative_start": pl.Int64,
    "relative_stop": pl.Int64,
    "pfam_ids": pl.List(pl.Utf8),
    "pfam_descriptions": pl.List(pl.Utf8),
    "interpro_ids": pl.List(pl.Utf8),
    "interpro_descriptions": pl.List(pl.Utf8),
    "sensory_domains": pl.List(pl.Utf8),
    "is_anchor": pl.Boolean,
    "is_regulator_candidate": pl.Boolean,
}


def _group_annotations(
    table: pl.DataFrame,
    *,
    id_column: str,
    description_column: str,
) -> pl.DataFrame:
    required = {"gene_accession", id_column, description_column}
    missing = required - set(table.columns)
    if missing:
        missing_text = ", ".join(sorted(missing))
        msg = f"annotation table is missing column(s): {missing_text}"
        raise ValueError(msg)

    if table.is_empty():
        return pl.DataFrame(
            schema={
                "gene_accession": pl.Utf8,
                f"{id_column}s": pl.List(pl.Utf8),
                f"{description_column}s": pl.List(pl.Utf8),
            },
        )

    return (
        table.select(
            pl.col("gene_accession").cast(pl.Utf8),
            pl.col(id_column).cast(pl.Utf8),
            pl.col(description_column).cast(pl.Utf8),
        )
        .group_by("gene_accession", maintain_order=True)
        .agg(
            pl.col(id_column)
            .drop_nulls()
            .unique(maintain_order=True)
            .alias(f"{id_column}s"),
            pl.col(description_column)
            .drop_nulls()
            .unique(maintain_order=True)
            .alias(f"{description_column}s"),
        )
    )


def annotate_domains(
    genes: pl.DataFrame,
    pfam_table: pl.DataFrame,
    interpro_table: pl.DataFrame,
) -> pl.DataFrame:
    """Populate Pfam and InterPro annotation list columns on a genes table."""

    pfam = _group_annotations(
        pfam_table,
        id_column="pfam_id",
        description_column="pfam_description",
    )
    interpro = _group_annotations(
        interpro_table,
        id_column="interpro_id",
        description_column="interpro_description",
    )

    annotated = genes.drop(
        [
            "pfam_ids",
            "pfam_descriptions",
            "interpro_ids",
            "interpro_descriptions",
        ],
    ).join(pfam, on="gene_accession", how="left").join(
        interpro,
        on="gene_accession",
        how="left",
    )
    empty_list = pl.lit([], dtype=pl.List(pl.Utf8))
    annotated = annotated.with_columns(
        pl.col("pfam_ids").fill_null(empty_list),
        pl.col("pfam_descriptions").fill_null(empty_list),
        pl.col("interpro_ids").fill_null(empty_list),
        pl.col("interpro_descriptions").fill_null(empty_list),
    )
    unannotated_by_locus = dict(
        annotated.filter(
            (pl.col("pfam_ids").list.len() == 0)
            & (pl.col("interpro_ids").list.len() == 0),
        )
        .group_by("locus_id", maintain_order=True)
        .len()
        .iter_rows(),
    )
    log.info("annotated domains", unannotated_by_locus=unannotated_by_locus)
    return validate(annotated, GenesSchema)
