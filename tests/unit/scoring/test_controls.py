from __future__ import annotations

from datetime import datetime

import polars as pl

from gasregnet.scoring.controls import sample_matched_controls


def _loci(ids: list[str], confidences: list[str]) -> pl.DataFrame:
    return pl.DataFrame(
        {
            "locus_id": ids,
            "analyte": ["CO"] * len(ids),
            "anchor_accession": [f"anchor_{idx}" for idx, _ in enumerate(ids)],
            "anchor_family": ["coxL"] * len(ids),
            "organism": ["OrgA"] * len(ids),
            "taxon_id": [1] * len(ids),
            "cluster_id": pl.Series([1] * len(ids), dtype=pl.Int32),
            "contig_id": ["contigA"] * len(ids),
            "window_size": pl.Series([10] * len(ids), dtype=pl.Int32),
            "is_boundary_truncated": [False] * len(ids),
            "marker_genes_present": pl.Series(
                [["coxL"]] * len(ids),
                dtype=pl.List(pl.Utf8),
            ),
            "accessory_genes_present": pl.Series(
                [[]] * len(ids),
                dtype=pl.List(pl.Utf8),
            ),
            "locus_score": [0.0] * len(ids),
            "locus_confidence": confidences,
            "taxonomic_context_score": [0.0] * len(ids),
            "operon_integrity_score": [0.0] * len(ids),
            "created_at": pl.Series(
                [datetime(2026, 4, 29)] * len(ids),
                dtype=pl.Datetime("us"),
            ),
        },
    )


def test_sample_matched_controls_is_deterministic() -> None:
    cases = _loci(["case1", "case2"], ["high", "high"])
    pool = _loci(
        ["case1", "case2", "ctrl1", "ctrl2", "ctrl3", "ctrl4", "high_ctrl"],
        ["high", "high", "low", "medium", "control", "low", "high"],
    )

    first = sample_matched_controls(cases, pool, ratio=(1, 2), seed=20260429)
    second = sample_matched_controls(cases, pool, ratio=(1, 2), seed=20260429)

    assert first["locus_id"].to_list() == second["locus_id"].to_list()
    assert first.height == 4
    assert "high_ctrl" not in first["locus_id"].to_list()
    assert set(first["locus_confidence"].to_list()) == {"control"}


def test_sample_matched_controls_returns_empty_when_no_match() -> None:
    cases = _loci(["case1"], ["high"]).with_columns(
        pl.lit("contigB").alias("contig_id"),
    )
    pool = _loci(["ctrl1"], ["low"])

    controls = sample_matched_controls(cases, pool, ratio=(1, 3), seed=1)

    assert controls.height == 0
