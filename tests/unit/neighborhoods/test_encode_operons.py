from __future__ import annotations

import polars as pl
import pytest

from gasregnet.neighborhoods.encode import encode_architectures
from gasregnet.neighborhoods.operons import (
    anchor_operon_integrity,
    infer_operon_membership,
)
from tests.unit.test_schemas import genes_frame, loci_frame


def _genes() -> pl.DataFrame:
    anchor = genes_frame().with_columns(
        pl.lit(101).alias("relative_start"),
        pl.lit(200).alias("relative_stop"),
    )
    upstream = genes_frame().with_columns(
        pl.lit("reg1").alias("gene_accession"),
        pl.lit(-1).cast(pl.Int32).alias("relative_index"),
        pl.lit(1).alias("relative_start"),
        pl.lit(80).alias("relative_stop"),
        pl.lit("+").alias("strand"),
        pl.lit("regulator").alias("functional_class"),
        pl.lit("one_component").alias("regulator_class"),
        pl.Series("sensory_domains", [["PAS"]], dtype=pl.List(pl.Utf8)),
        pl.lit(False).alias("is_anchor"),
        pl.lit(True).alias("is_regulator_candidate"),
    )
    downstream = genes_frame().with_columns(
        pl.lit("met1").alias("gene_accession"),
        pl.lit(1).cast(pl.Int32).alias("relative_index"),
        pl.lit(500).alias("relative_start"),
        pl.lit(600).alias("relative_stop"),
        pl.lit("-").alias("strand"),
        pl.lit("metabolic").alias("functional_class"),
        pl.lit(False).alias("is_anchor"),
        pl.lit(False).alias("is_regulator_candidate"),
    )
    return pl.concat([upstream, anchor, downstream])


def test_encode_architectures_centers_anchor_and_regulator_domains() -> None:
    architectures = encode_architectures(loci_frame(), _genes())

    assert architectures["architecture_string"].item() == (
        "[-1:one_component:PAS][0:coxL][+1:metabolic]"
    )


def test_infer_operon_membership_uses_distance_and_strand() -> None:
    memberships = infer_operon_membership(_genes(), max_intergenic_distance=25)

    upstream = memberships.filter(pl.col("gene_accession") == "reg1")
    anchor = memberships.filter(pl.col("gene_accession") == "anchor")
    downstream = memberships.filter(pl.col("gene_accession") == "met1")

    assert upstream["operon_id"].item() == anchor["operon_id"].item()
    assert upstream["in_anchor_operon"].item() is True
    assert downstream["in_anchor_operon"].item() is False


def test_anchor_operon_integrity_is_anchor_operon_fraction() -> None:
    memberships = infer_operon_membership(_genes(), max_intergenic_distance=25)

    integrity = anchor_operon_integrity(memberships)

    assert integrity["operon_integrity_score"].item() == pytest.approx(2 / 3)
