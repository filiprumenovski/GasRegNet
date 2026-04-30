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
) -> dict[str, tuple[list[str], list[str]]]:
    grouped: dict[str, tuple[list[str], list[str]]] = {}
    if table.is_empty():
        return grouped

    required = {"gene_accession", id_column, description_column}
    missing = required - set(table.columns)
    if missing:
        missing_text = ", ".join(sorted(missing))
        msg = f"annotation table is missing column(s): {missing_text}"
        raise ValueError(msg)

    for row in table.iter_rows(named=True):
        gene_accession = str(row["gene_accession"])
        ids, descriptions = grouped.setdefault(gene_accession, ([], []))
        annotation_id = row[id_column]
        description = row[description_column]
        if annotation_id is not None and str(annotation_id) not in ids:
            ids.append(str(annotation_id))
        if description is not None and str(description) not in descriptions:
            descriptions.append(str(description))
    return grouped


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

    rows: list[dict[str, object]] = []
    unannotated_by_locus: dict[str, int] = {}
    for row in genes.iter_rows(named=True):
        gene_accession = str(row["gene_accession"])
        pfam_ids, pfam_descriptions = pfam.get(gene_accession, ([], []))
        interpro_ids, interpro_descriptions = interpro.get(gene_accession, ([], []))
        updated = dict(row)
        updated["pfam_ids"] = pfam_ids
        updated["pfam_descriptions"] = pfam_descriptions
        updated["interpro_ids"] = interpro_ids
        updated["interpro_descriptions"] = interpro_descriptions
        if not pfam_ids and not interpro_ids:
            locus_id = str(row["locus_id"])
            unannotated_by_locus[locus_id] = unannotated_by_locus.get(locus_id, 0) + 1
        rows.append(updated)

    log.info("annotated domains", unannotated_by_locus=unannotated_by_locus)
    annotated = pl.DataFrame(
        rows,
        schema_overrides=GENES_SCHEMA_OVERRIDES,
    )
    return validate(annotated, GenesSchema)
