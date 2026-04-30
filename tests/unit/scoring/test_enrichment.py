from __future__ import annotations

import polars as pl
import pytest
from statsmodels.stats.multitest import multipletests  # type: ignore[import-untyped]

from gasregnet.scoring.enrichment import run_enrichment


def _genes(
    locus_ids: list[str],
    regulator_classes: list[str],
    sensory_domains: list[list[str]],
) -> pl.DataFrame:
    return pl.DataFrame(
        {
            "locus_id": locus_ids,
            "gene_accession": [f"{locus_id}_gene" for locus_id in locus_ids],
            "relative_index": pl.Series([-1] * len(locus_ids), dtype=pl.Int32),
            "relative_start": [0] * len(locus_ids),
            "relative_stop": [100] * len(locus_ids),
            "strand": ["+"] * len(locus_ids),
            "product_description": ["regulator"] * len(locus_ids),
            "pfam_ids": pl.Series(
                [["PF01047"]] * len(locus_ids),
                dtype=pl.List(pl.Utf8),
            ),
            "pfam_descriptions": pl.Series(
                [["domain"]] * len(locus_ids),
                dtype=pl.List(pl.Utf8),
            ),
            "interpro_ids": pl.Series([[]] * len(locus_ids), dtype=pl.List(pl.Utf8)),
            "interpro_descriptions": pl.Series(
                [[]] * len(locus_ids),
                dtype=pl.List(pl.Utf8),
            ),
            "functional_class": ["regulator"] * len(locus_ids),
            "regulator_class": regulator_classes,
            "sensory_domains": pl.Series(sensory_domains, dtype=pl.List(pl.Utf8)),
            "is_anchor": [False] * len(locus_ids),
            "is_regulator_candidate": [
                regulator_class != "none" for regulator_class in regulator_classes
            ],
        },
    )


def test_run_enrichment_recovers_known_odds_ratio_and_bh_q_values() -> None:
    case_genes = _genes(
        ["case1", "case2", "case3", "case4"],
        ["one_component", "one_component", "one_component", "none"],
        [["PAS"], ["PAS"], [], []],
    )
    control_genes = _genes(
        ["ctrl1", "ctrl2", "ctrl3", "ctrl4"],
        ["one_component", "none", "none", "none"],
        [[], [], [], []],
    )

    enrichment = run_enrichment(
        case_genes,
        control_genes,
        analyte="CO",
        case_definition="case",
        control_definition="control",
    )
    one_component = enrichment.filter(
        (pl.col("feature_type") == "regulator_class")
        & (pl.col("feature_name") == "one_component"),
    )

    assert one_component["n_case_with_feature"].item() == 3
    assert one_component["n_case_without_feature"].item() == 1
    assert one_component["n_control_with_feature"].item() == 1
    assert one_component["n_control_without_feature"].item() == 3
    assert one_component["odds_ratio"].item() == pytest.approx(9.0)

    expected_q = multipletests(
        enrichment["p_value"].to_list(),
        alpha=0.05,
        method="fdr_bh",
    )[1]
    assert enrichment["q_value"].to_list() == pytest.approx(expected_q.tolist())


def test_run_enrichment_returns_empty_schema_for_no_features() -> None:
    genes = _genes(["locus1"], ["none"], [[]])

    enrichment = run_enrichment(
        genes,
        genes,
        analyte="CO",
        case_definition="case",
        control_definition="control",
    )

    assert enrichment.height == 0
    assert set(enrichment.columns) == {
        "analyte",
        "feature_type",
        "feature_name",
        "case_definition",
        "control_definition",
        "n_case_with_feature",
        "n_case_without_feature",
        "n_control_with_feature",
        "n_control_without_feature",
        "odds_ratio",
        "p_value",
        "q_value",
        "interpretation",
    }
