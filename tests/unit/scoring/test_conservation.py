from __future__ import annotations

from datetime import datetime

import polars as pl

from gasregnet.archetypes.cluster import cluster_archetypes
from gasregnet.scoring.conservation import compute_conservation_scores


def _loci(prefix: str, genera: list[str], families: list[str]) -> pl.DataFrame:
    return pl.DataFrame(
        {
            "locus_id": [f"{prefix}_{index}" for index in range(len(genera))],
            "analyte": ["CO"] * len(genera),
            "anchor_accession": [f"anchor_{index}" for index in range(len(genera))],
            "anchor_family": ["coxL"] * len(genera),
            "organism": [f"{genus} species" for genus in genera],
            "taxon_id": list(range(100, 100 + len(genera))),
            "cluster_id": pl.Series([1] * len(genera), dtype=pl.Int32),
            "contig_id": [f"ctg_{index}" for index in range(len(genera))],
            "window_size": pl.Series([10] * len(genera), dtype=pl.Int32),
            "is_boundary_truncated": [False] * len(genera),
            "marker_genes_present": pl.Series(
                [["coxL"]] * len(genera),
                dtype=pl.List(pl.Utf8),
            ),
            "accessory_genes_present": pl.Series(
                [[]] * len(genera),
                dtype=pl.List(pl.Utf8),
            ),
            "locus_score": [6.0] * len(genera),
            "locus_confidence": ["high"] * len(genera),
            "taxonomic_context_score": [0.0] * len(genera),
            "operon_integrity_score": [0.0] * len(genera),
            "created_at": pl.Series(
                [datetime(2026, 4, 30)] * len(genera),
                dtype=pl.Datetime("us"),
            ),
            "genus": genera,
            "family": families,
            "phylum": ["Proteobacteria"] * len(genera),
        },
        schema_overrides={
            "marker_genes_present": pl.List(pl.Utf8),
            "accessory_genes_present": pl.List(pl.Utf8),
        },
    )


def _candidates(locus_ids: list[str]) -> pl.DataFrame:
    n = len(locus_ids)
    return pl.DataFrame(
        {
            "candidate_id": [f"{locus_id}::reg" for locus_id in locus_ids],
            "analyte": ["CO"] * n,
            "locus_id": locus_ids,
            "gene_accession": [f"reg_{index}" for index in range(n)],
            "organism": ["organism"] * n,
            "cluster_id": pl.Series([1] * n, dtype=pl.Int32),
            "relative_index": pl.Series([-1] * n, dtype=pl.Int32),
            "distance_nt": [100] * n,
            "position": ["upstream"] * n,
            "strand": ["+"] * n,
            "regulator_class": ["one_component"] * n,
            "dna_binding_domains": pl.Series([["PF01047"]] * n, dtype=pl.List(pl.Utf8)),
            "sensory_domains": pl.Series([["PAS"]] * n, dtype=pl.List(pl.Utf8)),
            "pfam_ids": pl.Series([["PF01047"]] * n, dtype=pl.List(pl.Utf8)),
            "interpro_ids": pl.Series([[]] * n, dtype=pl.List(pl.Utf8)),
            "archetype_id": [""] * n,
            "locus_score": [6.0] * n,
            "regulator_domain_score": [1.0] * n,
            "sensory_domain_score": [1.0] * n,
            "proximity_score": [0.5] * n,
            "archetype_conservation_score": [0.0] * n,
            "enrichment_score": [0.0] * n,
            "taxonomic_breadth_score": [0.0] * n,
            "structural_plausibility_score": [None] * n,
            "candidate_score": [8.0] * n,
            "candidate_score_q": [None] * n,
            "rationale": [""] * n,
        },
        schema_overrides={
            "structural_plausibility_score": pl.Float64,
            "candidate_score_q": pl.Float64,
        },
    )


def test_conservation_rewards_cross_genus_and_cross_family_architectures() -> None:
    broad_loci = _loci(
        "broad",
        ["GenusA", "GenusB", "GenusC", "GenusD", "GenusE", "GenusF", "GenusG"],
        ["FamA", "FamB", "FamC", "FamD", "FamE", "FamF", "FamG"],
    )
    narrow_loci = _loci(
        "narrow",
        ["GenusA"] * 7,
        ["FamA"] * 7,
    )
    broad_candidates = _candidates(broad_loci["locus_id"].to_list())
    narrow_candidates = _candidates(narrow_loci["locus_id"].to_list())

    broad = compute_conservation_scores(
        broad_candidates,
        cluster_archetypes(broad_loci, broad_candidates),
        broad_loci,
        min_loci_per_archetype=3,
    )
    narrow = compute_conservation_scores(
        narrow_candidates,
        cluster_archetypes(narrow_loci, narrow_candidates),
        narrow_loci,
        min_loci_per_archetype=3,
    )

    assert broad["archetype_conservation_score"].mean() > 0.9
    assert narrow["archetype_conservation_score"].mean() < 0.4
    assert broad["archetype_id"].str.len_chars().min() > 0
