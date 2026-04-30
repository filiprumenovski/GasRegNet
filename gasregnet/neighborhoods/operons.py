"""Operon inference helpers."""

from __future__ import annotations

from typing import Any

import polars as pl


def _intergenic_distance(left: dict[str, Any], right: dict[str, Any]) -> int:
    return max(0, int(right["relative_start"]) - int(left["relative_stop"]) - 1)


def _same_operon(
    left: dict[str, Any],
    right: dict[str, Any],
    max_intergenic_distance: int,
) -> bool:
    return (
        str(left["strand"]) == str(right["strand"])
        and _intergenic_distance(left, right) <= max_intergenic_distance
    )


def infer_operon_membership(
    genes: pl.DataFrame,
    *,
    max_intergenic_distance: int = 200,
) -> pl.DataFrame:
    """Assign deterministic operon groups within each locus."""

    rows: list[dict[str, Any]] = []
    for (locus_id,), locus_genes in genes.group_by("locus_id", maintain_order=True):
        ordered = list(
            locus_genes.sort(["relative_start", "relative_stop"]).iter_rows(named=True),
        )
        group = 0
        group_by_accession: dict[str, int] = {}
        previous: dict[str, Any] | None = None
        for gene in ordered:
            if previous is not None and not _same_operon(
                previous,
                gene,
                max_intergenic_distance,
            ):
                group += 1
            accession = str(gene["gene_accession"])
            group_by_accession[accession] = group
            previous = gene

        anchor_groups = {
            group_by_accession[str(gene["gene_accession"])]
            for gene in ordered
            if int(gene["relative_index"]) == 0
        }
        for gene in ordered:
            accession = str(gene["gene_accession"])
            operon_group = group_by_accession[accession]
            row = dict(gene)
            row["operon_id"] = f"{locus_id}:operon_{operon_group + 1:02d}"
            row["in_anchor_operon"] = operon_group in anchor_groups
            rows.append(row)

    return pl.DataFrame(rows) if rows else genes.with_columns(
        pl.lit("").alias("operon_id"),
        pl.lit(False).alias("in_anchor_operon"),
    )


def anchor_operon_integrity(genes: pl.DataFrame) -> pl.DataFrame:
    """Summarize the fraction of genes assigned to each anchor operon."""

    if "in_anchor_operon" not in genes.columns:
        genes = infer_operon_membership(genes)
    return genes.group_by("locus_id").agg(
        pl.col("in_anchor_operon").mean().alias("operon_integrity_score"),
    )
