from __future__ import annotations

from datetime import datetime

import polars as pl
import pytest

from gasregnet.archetypes.cluster import architecture_distance, cluster_archetypes


def _loci() -> pl.DataFrame:
    return pl.DataFrame(
        {
            "locus_id": ["locus_a", "locus_b", "locus_c"],
            "analyte": ["CO", "CO", "CO"],
            "anchor_accession": ["a", "b", "c"],
            "anchor_family": ["coxL", "coxL", "coxL"],
            "organism": ["Org1", "Org2", "Org3"],
            "taxon_id": [1, 2, 2],
            "cluster_id": pl.Series([1, 1, 1], dtype=pl.Int32),
            "contig_id": ["ctg1", "ctg2", "ctg3"],
            "window_size": pl.Series([10, 10, 10], dtype=pl.Int32),
            "is_boundary_truncated": [False, False, False],
            "marker_genes_present": pl.Series(
                [["coxL"], ["coxL"], ["coxL"]],
                dtype=pl.List(pl.Utf8),
            ),
            "accessory_genes_present": pl.Series(
                [[], [], []],
                dtype=pl.List(pl.Utf8),
            ),
            "locus_score": [6.0, 5.0, 4.0],
            "locus_confidence": ["high", "high", "medium"],
            "taxonomic_context_score": [0.0, 0.0, 0.0],
            "operon_integrity_score": [0.0, 0.0, 0.0],
            "created_at": pl.Series(
                [datetime(2026, 4, 29)] * 3,
                dtype=pl.Datetime("us"),
            ),
        },
    )


def _candidates() -> pl.DataFrame:
    return pl.DataFrame(
        {
            "candidate_id": ["ca", "cb", "cc"],
            "analyte": ["CO", "CO", "CO"],
            "locus_id": ["locus_a", "locus_b", "locus_c"],
            "gene_accession": ["rega", "regb", "regc"],
            "organism": ["Org1", "Org2", "Org3"],
            "cluster_id": pl.Series([1, 1, 1], dtype=pl.Int32),
            "relative_index": pl.Series([-1, -1, -2], dtype=pl.Int32),
            "distance_nt": [100, 120, 300],
            "position": ["upstream", "upstream", "upstream"],
            "strand": ["+", "+", "+"],
            "regulator_class": ["one_component", "one_component", "two_component_rr"],
            "dna_binding_domains": pl.Series(
                [["PF01047"], ["PF01047"], ["PF00072"]],
                dtype=pl.List(pl.Utf8),
            ),
            "sensory_domains": pl.Series(
                [["PAS"], ["PAS"], []],
                dtype=pl.List(pl.Utf8),
            ),
            "pfam_ids": pl.Series(
                [["PF01047"], ["PF01047"], ["PF00072"]],
                dtype=pl.List(pl.Utf8),
            ),
            "interpro_ids": pl.Series([[], [], []], dtype=pl.List(pl.Utf8)),
            "archetype_id": ["", "", ""],
            "locus_score": [6.0, 5.0, 4.0],
            "regulator_domain_score": [1.0, 1.0, 1.0],
            "sensory_domain_score": [1.0, 1.0, 0.0],
            "proximity_score": [0.5, 0.5, 0.333],
            "archetype_conservation_score": [0.0, 0.0, 0.0],
            "enrichment_score": [0.0, 0.0, 0.0],
            "taxonomic_breadth_score": [0.0, 0.0, 0.0],
            "structural_plausibility_score": [None, None, None],
            "candidate_score": [8.0, 7.0, 5.0],
            "candidate_score_q": [None, None, None],
            "rationale": ["a", "b", "c"],
        },
        schema_overrides={
            "structural_plausibility_score": pl.Float64,
            "candidate_score_q": pl.Float64,
        },
    )


def test_cluster_archetypes_groups_identical_architectures() -> None:
    archetypes = cluster_archetypes(_loci(), _candidates())

    assert archetypes.height == 2
    first = archetypes.filter(pl.col("archetype_id") == "arch_0001")
    assert first["n_loci"].item() == 2
    assert first["n_taxa"].item() == 2
    assert first["representative_locus_id"].item() == "locus_a"
    assert first["dominant_regulator_class"].item() == "one_component"
    assert first["mean_candidate_score"].item() == pytest.approx(7.5)


def test_cluster_archetypes_is_deterministic() -> None:
    first = cluster_archetypes(_loci(), _candidates())
    second = cluster_archetypes(_loci(), _candidates())

    assert first.equals(second)


def test_architecture_distance_weights_near_anchor_more_heavily() -> None:
    near_change = architecture_distance(
        "[-1:one_component:PAS][0:coxL]",
        "[-1:two_component_rr:none][0:coxL]",
    )
    far_change = architecture_distance(
        "[-5:one_component:PAS][0:coxL]",
        "[-5:two_component_rr:none][0:coxL]",
    )

    assert near_change > far_change


def test_cluster_archetypes_threshold_can_merge_related_architectures() -> None:
    archetypes = cluster_archetypes(_loci(), _candidates(), distance_threshold=1.0)

    assert archetypes.height == 1
    assert archetypes["n_loci"].item() == 3
