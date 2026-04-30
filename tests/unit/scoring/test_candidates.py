from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest

from gasregnet.config import load_config
from gasregnet.scoring.candidates import score_candidates
from tests.unit.test_schemas import loci_frame


def _genes() -> pl.DataFrame:
    return pl.DataFrame(
        {
            "locus_id": ["CO_1_anchor", "CO_1_anchor", "CO_1_anchor"],
            "gene_accession": ["near_pas", "far_plain", "medium_two_domain"],
            "relative_index": pl.Series([-1, 4, 2], dtype=pl.Int32),
            "relative_start": [-150, 1000, 500],
            "relative_stop": [-50, 1200, 650],
            "strand": ["+", "+", "-"],
            "product_description": ["regulator", "regulator", "regulator"],
            "pfam_ids": pl.Series(
                [["PF01047", "PF00989"], ["PF01047"], ["PF01047", "PF00989"]],
                dtype=pl.List(pl.Utf8),
            ),
            "pfam_descriptions": pl.Series(
                [["MarR", "PAS"], ["MarR"], ["MarR", "PAS"]],
                dtype=pl.List(pl.Utf8),
            ),
            "interpro_ids": pl.Series(
                [["IPR1"], ["IPR2"], ["IPR3"]],
                dtype=pl.List(pl.Utf8),
            ),
            "interpro_descriptions": pl.Series(
                [["domain"], ["domain"], ["domain"]],
                dtype=pl.List(pl.Utf8),
            ),
            "functional_class": ["regulator", "regulator", "regulator"],
            "regulator_class": ["one_component", "one_component", "one_component"],
            "sensory_domains": pl.Series(
                [["PAS"], [], ["PAS", "GAF"]],
                dtype=pl.List(pl.Utf8),
            ),
            "is_anchor": [False, False, False],
            "is_regulator_candidate": [True, True, True],
        },
    )


def _enrichment() -> pl.DataFrame:
    return pl.DataFrame(
        {
            "analyte": ["CO", "CO"],
            "feature_type": ["regulator_class", "sensory_domain"],
            "feature_name": ["one_component", "PAS"],
            "case_definition": ["case", "case"],
            "control_definition": ["control", "control"],
            "n_case_with_feature": [3, 2],
            "n_case_without_feature": [0, 1],
            "n_control_with_feature": [1, 0],
            "n_control_without_feature": [2, 3],
            "odds_ratio": [4.0, 8.0],
            "p_value": [0.05, 0.01],
            "q_value": [0.04, 0.02],
            "interpretation": ["enriched", "enriched"],
        },
    )


def test_score_candidates_reproduces_hand_computed_ranking() -> None:
    config = load_config(Path("configs"))
    loci = loci_frame().with_columns(pl.lit(6.0).alias("locus_score"))

    candidates = score_candidates(loci, _genes(), config.scoring, _enrichment())
    ranking = candidates["gene_accession"].to_list()

    assert ranking == ["medium_two_domain", "near_pas", "far_plain"]

    top = candidates.filter(pl.col("gene_accession") == "medium_two_domain")
    expected = 6.0 + 1.5 + 4.0 + (1.0 / 3.0) + 6.0
    assert top["candidate_score"].item() == pytest.approx(expected)
    assert top["regulator_domain_score"].item() == 1.0
    assert top["sensory_domain_score"].item() == 2.0
    assert top["enrichment_score"].item() == pytest.approx(3.0)
    assert top["candidate_score_q"].item() == 0.02


def test_score_candidates_populates_required_fields_and_rationale() -> None:
    config = load_config(Path("configs"))
    candidates = score_candidates(loci_frame(), _genes(), config.scoring)

    row = candidates.filter(pl.col("gene_accession") == "near_pas")

    assert row["candidate_id"].item() == "CO_1_anchor::near_pas"
    assert row["position"].item() == "upstream"
    assert row["distance_nt"].item() == 50
    assert row["dna_binding_domains"].to_list()[0] == ["PF01047"]
    assert row["rationale"].item()


def test_score_candidates_returns_empty_schema_when_no_candidates() -> None:
    config = load_config(Path("configs"))
    genes = _genes().with_columns(pl.lit(False).alias("is_regulator_candidate"))

    candidates = score_candidates(loci_frame(), genes, config.scoring)

    assert candidates.height == 0
    assert "candidate_score" in candidates.columns
