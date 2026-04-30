from __future__ import annotations

from datetime import datetime

import polars as pl

from gasregnet.annotation.roles import (
    assign_sensor_roles,
    build_sensor_regulator_pairs,
)
from gasregnet.config import load_config


def _genes(pfams: list[list[str]]) -> pl.DataFrame:
    return pl.DataFrame(
        {
            "locus_id": ["locus1"] * len(pfams),
            "gene_accession": [f"gene{i}" for i in range(len(pfams))],
            "relative_index": pl.Series(range(1, len(pfams) + 1), dtype=pl.Int32),
            "relative_start": [i * 100 for i in range(len(pfams))],
            "relative_stop": [i * 100 + 80 for i in range(len(pfams))],
            "strand": ["+"] * len(pfams),
            "product_description": ["annotated"] * len(pfams),
            "pfam_ids": pl.Series(pfams, dtype=pl.List(pl.Utf8)),
            "pfam_descriptions": pl.Series([[]] * len(pfams), dtype=pl.List(pl.Utf8)),
            "interpro_ids": pl.Series([[]] * len(pfams), dtype=pl.List(pl.Utf8)),
            "interpro_descriptions": pl.Series(
                [[]] * len(pfams),
                dtype=pl.List(pl.Utf8),
            ),
            "functional_class": ["unknown"] * len(pfams),
            "regulator_class": ["none"] * len(pfams),
            "sensory_domains": pl.Series([[]] * len(pfams), dtype=pl.List(pl.Utf8)),
            "is_anchor": [False] * len(pfams),
            "is_regulator_candidate": [False] * len(pfams),
        },
    )


def _loci() -> pl.DataFrame:
    return pl.DataFrame(
        {
            "locus_id": ["locus1"],
            "analyte": ["CO"],
            "anchor_accession": ["anchor"],
            "anchor_family": ["fix"],
            "organism": ["Bradyrhizobium diazoefficiens"],
            "taxon_id": [224911],
            "cluster_id": pl.Series([1], dtype=pl.Int32),
            "contig_id": ["ctg"],
            "window_size": pl.Series([3], dtype=pl.Int32),
            "is_boundary_truncated": [False],
            "marker_genes_present": [["fix"]],
            "accessory_genes_present": [[]],
            "locus_score": [0.0],
            "locus_confidence": ["low"],
            "taxonomic_context_score": [0.0],
            "operon_integrity_score": [0.0],
            "created_at": pl.Series([datetime(2026, 4, 30)], dtype=pl.Datetime("us")),
        },
        schema_overrides={
            "marker_genes_present": pl.List(pl.Utf8),
            "accessory_genes_present": pl.List(pl.Utf8),
        },
    )


def _assign(pfams: list[list[str]]) -> pl.DataFrame:
    config = load_config("configs")
    return assign_sensor_roles(
        _genes(pfams),
        regulator_families=config.regulator_families,
        sensory_domain_catalog=config.sensory_domains,
        paired_evidence_rules=config.paired_evidence,
    )


def test_assign_sensor_roles_classifies_cooa_as_both() -> None:
    assigned = _assign([["PF00027", "PF00325"]])

    assert assigned["regulator_class"].item() == "one_component"
    assert assigned["sensor_role"].item() == "both"
    assert assigned["is_regulator_candidate"].item() is True


def test_assign_sensor_roles_classifies_fixl_as_sensor_and_fixj_as_regulator() -> None:
    assigned = _assign([["PF00512", "PF00989"], ["PF00072"]])

    assert assigned["sensor_role"].to_list() == ["sensor", "regulator"]
    assert assigned["is_regulator_candidate"].to_list() == [False, True]


def test_build_sensor_regulator_pairs_links_fixl_fixj() -> None:
    config = load_config("configs")
    assigned = assign_sensor_roles(
        _genes([["PF00512", "PF00989"], ["PF00072"]]),
        regulator_families=config.regulator_families,
        sensory_domain_catalog=config.sensory_domains,
        paired_evidence_rules=config.paired_evidence,
    )

    pairs = build_sensor_regulator_pairs(
        assigned,
        _loci(),
        sensory_domain_catalog=config.sensory_domains,
        paired_evidence_rules=config.paired_evidence,
    )

    assert pairs.height == 1
    assert pairs["hk_gene_accession"].item() == "gene0"
    assert pairs["rr_gene_accession"].item() == "gene1"
    assert pairs["hk_sensor_domains"].to_list()[0] == ["PAS_generic"]


def test_assign_sensor_roles_classifies_rcom_and_arsr_as_both() -> None:
    assigned = _assign([["PF00027", "PF00325", "PF00989"], ["PF01022"]])

    assert assigned["sensor_role"].to_list() == ["both", "both"]
