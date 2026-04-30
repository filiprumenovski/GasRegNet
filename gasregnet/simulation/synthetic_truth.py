"""Synthetic-truth genomes for calibration testing."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from random import Random

import polars as pl

from gasregnet.schemas import GenesSchema, LociSchema, validate

COOA_OPERATOR = "TTGTCA-N6-TGACAA"


@dataclass(frozen=True)
class SyntheticTruthCorpus:
    """Known-truth frames for end-to-end calibration tests."""

    loci: pl.DataFrame
    genes: pl.DataFrame
    ground_truth: pl.DataFrame


def _regulator_class(index: int, has_anchor: bool, rng: Random, noise: float) -> str:
    if has_anchor:
        return "one_component" if rng.random() >= noise else "none"
    return "one_component" if rng.random() < noise else "none"


def simulate_synthetic_truth_corpus(
    *,
    n_genomes: int = 24,
    positive_fraction: float = 0.5,
    phylogenetic_clumping: float = 0.7,
    annotation_noise: float = 0.1,
    seed: int = 20260430,
) -> SyntheticTruthCorpus:
    """Generate synthetic genomes with known planted regulator-anchor truth.

    Positive genomes carry a CooA-like regulator upstream of a synthetic CO anchor
    and a planted operator sequence. Negative genomes carry a decoy locus without
    the planted regulator-anchor relation. ``phylogenetic_clumping`` controls how
    strongly positives are concentrated in adjacent synthetic clades.
    """

    if not 0.0 <= positive_fraction <= 1.0:
        raise ValueError("positive_fraction must be in [0, 1]")
    if not 0.0 <= phylogenetic_clumping <= 1.0:
        raise ValueError("phylogenetic_clumping must be in [0, 1]")
    if not 0.0 <= annotation_noise <= 1.0:
        raise ValueError("annotation_noise must be in [0, 1]")

    rng = Random(seed)
    n_positive = round(n_genomes * positive_fraction)
    clade_size = max(1, round(max(1, n_genomes // 4) * phylogenetic_clumping))
    positive_indexes = set(range(min(n_positive, clade_size)))
    while len(positive_indexes) < n_positive:
        positive_indexes.add(rng.randrange(n_genomes))

    loci_rows: list[dict[str, object]] = []
    gene_rows: list[dict[str, object]] = []
    truth_rows: list[dict[str, object]] = []
    created_at = datetime(2026, 4, 30)

    for index in range(n_genomes):
        has_truth = index in positive_indexes
        analyte = "CO" if has_truth else "CN"
        anchor_family = "synthetic_coxL" if has_truth else "synthetic_cydA"
        organism = f"Synthetic bacterium {index:03d}"
        locus_id = f"SYN_{index:03d}_{analyte.lower()}_anchor"
        cluster_id = index // max(1, n_genomes // 4)
        loci_rows.append(
            {
                "locus_id": locus_id,
                "analyte": analyte,
                "anchor_accession": f"SYNANCHOR_{index:03d}",
                "anchor_family": anchor_family,
                "organism": organism,
                "taxon_id": 9_000_000 + index,
                "cluster_id": cluster_id,
                "contig_id": f"synthetic_contig_{index:03d}",
                "window_size": 8,
                "is_boundary_truncated": False,
                "marker_genes_present": [anchor_family],
                "accessory_genes_present": ["synthetic_coxM"] if has_truth else [],
                "locus_score": 0.0,
                "locus_confidence": "low",
                "taxonomic_context_score": 0.0,
                "operon_integrity_score": 0.0,
                "created_at": created_at,
                "phylum": f"synthetic_phylum_{cluster_id}",
                "class": f"synthetic_class_{cluster_id}",
                "order": f"synthetic_order_{cluster_id}",
                "family": f"synthetic_family_{cluster_id}",
                "genus": f"Synthgenus{cluster_id}",
            },
        )
        regulator_class = _regulator_class(index, has_truth, rng, annotation_noise)
        sensory_domains = ["PAS", "CooA_heme"] if regulator_class != "none" else []
        pfam_ids = ["PF00196", "PF00989"] if regulator_class != "none" else []
        gene_rows.extend(
            [
                {
                    "locus_id": locus_id,
                    "gene_accession": f"SYNREG_{index:03d}",
                    "relative_index": -1,
                    "relative_start": -420,
                    "relative_stop": -90,
                    "strand": "+",
                    "product_description": (
                        "synthetic CooA-like regulator with planted operator "
                        f"{COOA_OPERATOR}"
                        if has_truth
                        else "synthetic decoy regulator"
                    ),
                    "pfam_ids": pfam_ids,
                    "pfam_descriptions": ["HTH", "PAS"] if pfam_ids else [],
                    "interpro_ids": ["IPR000000"] if pfam_ids else [],
                    "interpro_descriptions": ["synthetic regulator"]
                    if pfam_ids
                    else [],
                    "functional_class": "regulator"
                    if regulator_class != "none"
                    else "unknown",
                    "regulator_class": regulator_class,
                    "sensory_domains": sensory_domains,
                    "is_anchor": False,
                    "is_regulator_candidate": regulator_class != "none",
                },
                {
                    "locus_id": locus_id,
                    "gene_accession": f"SYNANCHOR_{index:03d}",
                    "relative_index": 0,
                    "relative_start": 1,
                    "relative_stop": 1200,
                    "strand": "+",
                    "product_description": (
                        "synthetic carbon monoxide dehydrogenase"
                        if has_truth
                        else "synthetic cyanide oxidase"
                    ),
                    "pfam_ids": ["PF02738"] if has_truth else ["PF01654"],
                    "pfam_descriptions": ["coxL"] if has_truth else ["cydA"],
                    "interpro_ids": ["IPR000000"],
                    "interpro_descriptions": ["synthetic anchor"],
                    "functional_class": "anchor",
                    "regulator_class": "none",
                    "sensory_domains": [],
                    "is_anchor": True,
                    "is_regulator_candidate": False,
                },
            ],
        )
        truth_rows.append(
            {
                "locus_id": locus_id,
                "gene_accession": f"SYNREG_{index:03d}",
                "organism": organism,
                "analyte": analyte,
                "is_true_regulator": has_truth,
                "operator_motif": COOA_OPERATOR if has_truth else "",
                "phylogenetic_clump": cluster_id,
                "annotation_noisy": regulator_class == "none" and has_truth,
            },
        )

    loci = validate(
        pl.DataFrame(
            loci_rows,
            schema_overrides={
                "cluster_id": pl.Int32,
                "window_size": pl.Int32,
                "created_at": pl.Datetime("us"),
            },
        ),
        LociSchema,
    )
    genes = validate(
        pl.DataFrame(
            gene_rows,
            schema_overrides={
                "relative_index": pl.Int32,
                "pfam_ids": pl.List(pl.Utf8),
                "pfam_descriptions": pl.List(pl.Utf8),
                "interpro_ids": pl.List(pl.Utf8),
                "interpro_descriptions": pl.List(pl.Utf8),
                "sensory_domains": pl.List(pl.Utf8),
            },
        ),
        GenesSchema,
    )
    ground_truth = pl.DataFrame(truth_rows)
    return SyntheticTruthCorpus(loci=loci, genes=genes, ground_truth=ground_truth)
