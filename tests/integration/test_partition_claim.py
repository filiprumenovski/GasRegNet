from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import cast

import polars as pl

from gasregnet.config import load_config
from gasregnet.scoring.candidates import score_candidates
from gasregnet.scoring.enrichment import run_enrichment
from gasregnet.scoring.partition import (
    chemistry_partition_outcome,
    family_chemistry_table,
    write_partition_outcome,
)


def _loci(analyte: str, locus_ids: list[str]) -> pl.DataFrame:
    anchor_family = "coxL" if analyte == "CO" else "cydA"
    return pl.DataFrame(
        {
            "locus_id": locus_ids,
            "analyte": [analyte] * len(locus_ids),
            "anchor_accession": [f"{locus_id}_anchor" for locus_id in locus_ids],
            "anchor_family": [anchor_family] * len(locus_ids),
            "organism": [f"{analyte}_Org"] * len(locus_ids),
            "taxon_id": [1 if analyte == "CO" else 2] * len(locus_ids),
            "cluster_id": pl.Series([1] * len(locus_ids), dtype=pl.Int32),
            "contig_id": [f"{analyte}_contig"] * len(locus_ids),
            "window_size": pl.Series([10] * len(locus_ids), dtype=pl.Int32),
            "is_boundary_truncated": [False] * len(locus_ids),
            "marker_genes_present": pl.Series(
                [[anchor_family]] * len(locus_ids),
                dtype=pl.List(pl.Utf8),
            ),
            "accessory_genes_present": pl.Series(
                [[]] * len(locus_ids),
                dtype=pl.List(pl.Utf8),
            ),
            "locus_score": [6.0] * len(locus_ids),
            "locus_confidence": ["high"] * len(locus_ids),
            "taxonomic_context_score": [0.0] * len(locus_ids),
            "operon_integrity_score": [0.0] * len(locus_ids),
            "created_at": pl.Series(
                [datetime(2026, 4, 29)] * len(locus_ids),
                dtype=pl.Datetime("us"),
            ),
        },
    )


def _genes(
    locus_ids: list[str],
    regulator_class: str,
    sensory_domains: list[list[str]],
) -> pl.DataFrame:
    return pl.DataFrame(
        {
            "locus_id": locus_ids,
            "gene_accession": [f"{locus_id}_reg" for locus_id in locus_ids],
            "relative_index": pl.Series([-1] * len(locus_ids), dtype=pl.Int32),
            "relative_start": [-200] * len(locus_ids),
            "relative_stop": [-50] * len(locus_ids),
            "strand": ["+"] * len(locus_ids),
            "product_description": ["regulator"] * len(locus_ids),
            "pfam_ids": pl.Series(
                [["PF01047"]] * len(locus_ids),
                dtype=pl.List(pl.Utf8),
            ),
            "pfam_descriptions": pl.Series(
                [["regulator"]] * len(locus_ids),
                dtype=pl.List(pl.Utf8),
            ),
            "interpro_ids": pl.Series([[]] * len(locus_ids), dtype=pl.List(pl.Utf8)),
            "interpro_descriptions": pl.Series(
                [[]] * len(locus_ids),
                dtype=pl.List(pl.Utf8),
            ),
            "functional_class": ["regulator"] * len(locus_ids),
            "regulator_class": [regulator_class] * len(locus_ids),
            "sensory_domains": pl.Series(sensory_domains, dtype=pl.List(pl.Utf8)),
            "is_anchor": [False] * len(locus_ids),
            "is_regulator_candidate": [True] * len(locus_ids),
        },
    )


def test_partition_claim_writes_outcome_json() -> None:
    config = load_config(Path("configs"))
    co_loci = _loci("CO", ["co1", "co2", "co3", "co4"])
    cn_loci = _loci("CN", ["cn1", "cn2", "cn3", "cn4"])
    co_case = _genes(co_loci["locus_id"].to_list(), "one_component", [["PAS"]] * 4)
    cn_case = _genes(cn_loci["locus_id"].to_list(), "two_component_rr", [["GAF"]] * 4)
    co_control = _genes(["co_ctrl1", "co_ctrl2"], "none", [[], []])
    cn_control = _genes(["cn_ctrl1", "cn_ctrl2"], "none", [[], []])

    co_enrichment = run_enrichment(
        co_case,
        co_control,
        analyte="CO",
        case_definition="CO synthetic case",
        control_definition="CO synthetic control",
    )
    cn_enrichment = run_enrichment(
        cn_case,
        cn_control,
        analyte="CN",
        case_definition="CN synthetic case",
        control_definition="CN synthetic control",
    )
    assert co_enrichment.height > 0
    assert cn_enrichment.height > 0

    candidates = pl.concat(
        [
            score_candidates(co_loci, co_case, config.scoring, co_enrichment),
            score_candidates(cn_loci, cn_case, config.scoring, cn_enrichment),
        ],
    )
    outcome = chemistry_partition_outcome(candidates)
    out_path = write_partition_outcome(
        outcome,
        Path("results/headline/partition_outcome.json"),
    )

    assert out_path.exists()
    table = family_chemistry_table(candidates)
    assert table["n_candidates"].sum() == candidates.height
    assert (
        outcome["reason"]
        == "chi-square test on analyte-by-primary-sensory-chemistry counts"
    )
    assert isinstance(outcome["partition_holds"], bool)
    assert float(cast(float, outcome["p_value"])) >= 0.0
