from __future__ import annotations

from pathlib import Path

import polars as pl

from gasregnet.config import load_config
from gasregnet.scoring.loci import score_loci
from tests.unit.test_schemas import loci_frame


def test_score_loci_computes_weighted_total() -> None:
    config = load_config(Path("configs"))
    loci = loci_frame().with_columns(
        pl.Series("marker_genes_present", [["coxL"]], dtype=pl.List(pl.Utf8)),
        pl.Series(
            "accessory_genes_present",
            [["coxM", "coxS"]],
            dtype=pl.List(pl.Utf8),
        ),
        pl.lit(1.0).alias("operon_integrity_score"),
        pl.lit(0.5).alias("taxonomic_context_score"),
        pl.lit(False).alias("is_boundary_truncated"),
    )

    scored = score_loci(loci, config.scoring)

    expected = 3.0 + 2.0 + 1.5 + 1.0 + 0.25 + 0.5
    assert scored["locus_score"].item() == expected
    assert scored["anchor_marker_score"].item() == 1.0
    assert scored["accessory_marker_score"].item() == 2.0
    assert scored["locus_confidence"].item() == "high"


def test_score_loci_is_deterministic() -> None:
    config = load_config(Path("configs"))

    first = score_loci(loci_frame(), config.scoring)
    second = score_loci(loci_frame(), config.scoring)

    assert first.select(pl.exclude("created_at")).equals(
        second.select(pl.exclude("created_at")),
    )


def test_score_loci_assigns_confidence_thresholds() -> None:
    config = load_config(Path("configs"))
    loci = pl.concat(
        [
            loci_frame().with_columns(
                pl.lit("high_case").alias("locus_id"),
                pl.Series("marker_genes_present", [["coxL"]], dtype=pl.List(pl.Utf8)),
                pl.Series(
                    "accessory_genes_present",
                    [["coxM"]],
                    dtype=pl.List(pl.Utf8),
                ),
                pl.lit(1.0).alias("operon_integrity_score"),
                pl.lit(0.0).alias("taxonomic_context_score"),
                pl.lit(False).alias("is_boundary_truncated"),
            ),
            loci_frame().with_columns(
                pl.lit("medium_case").alias("locus_id"),
                pl.Series("marker_genes_present", [["coxL"]], dtype=pl.List(pl.Utf8)),
                pl.Series("accessory_genes_present", [[]], dtype=pl.List(pl.Utf8)),
                pl.lit(0.0).alias("operon_integrity_score"),
                pl.lit(0.0).alias("taxonomic_context_score"),
                pl.lit(False).alias("is_boundary_truncated"),
            ),
            loci_frame().with_columns(
                pl.lit("control_case").alias("locus_id"),
                pl.Series("marker_genes_present", [[]], dtype=pl.List(pl.Utf8)),
                pl.Series("accessory_genes_present", [[]], dtype=pl.List(pl.Utf8)),
                pl.lit("").alias("anchor_accession"),
                pl.lit(0.0).alias("operon_integrity_score"),
                pl.lit(0.0).alias("taxonomic_context_score"),
                pl.lit(True).alias("is_boundary_truncated"),
            ),
        ],
    )

    scored = score_loci(loci, config.scoring)
    confidence = dict(zip(scored["locus_id"], scored["locus_confidence"], strict=True))

    assert confidence == {
        "high_case": "high",
        "medium_case": "medium",
        "control_case": "control",
    }


def test_score_loci_penalizes_boundary_truncation() -> None:
    config = load_config(Path("configs"))
    complete = loci_frame().with_columns(pl.lit(False).alias("is_boundary_truncated"))
    truncated = loci_frame().with_columns(pl.lit(True).alias("is_boundary_truncated"))

    complete_score = score_loci(complete, config.scoring)["locus_score"].item()
    truncated_score = score_loci(truncated, config.scoring)["locus_score"].item()

    assert complete_score - truncated_score == 0.5
