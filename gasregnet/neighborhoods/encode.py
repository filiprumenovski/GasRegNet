"""Neighborhood architecture encoding."""

from __future__ import annotations

from typing import Any

import polars as pl


def _clean_label(value: object) -> str:
    label = str(value).strip().replace(" ", "_")
    return label or "unknown"


def _gene_token(row: dict[str, Any], anchor_family: str) -> str:
    relative_index = int(row["relative_index"])
    if relative_index == 0:
        label = anchor_family
    elif row.get("is_regulator_candidate") is True:
        sensory_domains = "-".join(row.get("sensory_domains") or []) or "none"
        regulator_class = _clean_label(row.get("regulator_class", "regulator"))
        label = f"{regulator_class}:{sensory_domains}"
    else:
        label = _clean_label(row.get("functional_class", "unknown"))
    sign = "+" if relative_index > 0 else ""
    return f"[{sign}{relative_index}:{label}]"


def encode_locus_architecture(
    locus: dict[str, Any],
    genes: pl.DataFrame,
) -> str:
    """Encode one centered neighborhood architecture string."""

    locus_id = str(locus["locus_id"])
    anchor_family = str(locus["anchor_family"])
    locus_genes = genes.filter(pl.col("locus_id") == locus_id).sort("relative_index")
    tokens = [
        _gene_token(row, anchor_family)
        for row in locus_genes.iter_rows(named=True)
    ]
    return "".join(tokens)


def encode_architectures(loci: pl.DataFrame, genes: pl.DataFrame) -> pl.DataFrame:
    """Return one architecture string per locus."""

    rows = [
        {
            "locus_id": str(locus["locus_id"]),
            "architecture_string": encode_locus_architecture(locus, genes),
        }
        for locus in loci.iter_rows(named=True)
    ]
    return pl.DataFrame(
        rows,
        schema_overrides={
            "locus_id": pl.Utf8,
            "architecture_string": pl.Utf8,
        },
    )
