"""Candidate regulator scoring."""

from __future__ import annotations

import math
from typing import Any, cast

import polars as pl

from gasregnet.config import ScoringConfig
from gasregnet.schemas import RegulatorCandidatesSchema, validate

DNA_BINDING_PFAMS = {
    "PF00126",
    "PF00196",
    "PF00325",
    "PF00356",
    "PF00376",
    "PF00392",
    "PF00440",
    "PF01022",
    "PF01047",
    "PF01381",
    "PF03466",
    "PF04542",
    "PF04545",
    "PF12833",
}

CANDIDATE_SCHEMA: dict[str, Any] = {
    "candidate_id": pl.Utf8,
    "analyte": pl.Utf8,
    "locus_id": pl.Utf8,
    "gene_accession": pl.Utf8,
    "organism": pl.Utf8,
    "cluster_id": pl.Int32,
    "relative_index": pl.Int32,
    "distance_nt": pl.Int64,
    "position": pl.Utf8,
    "strand": pl.Utf8,
    "regulator_class": pl.Utf8,
    "dna_binding_domains": pl.List(pl.Utf8),
    "sensory_domains": pl.List(pl.Utf8),
    "pfam_ids": pl.List(pl.Utf8),
    "interpro_ids": pl.List(pl.Utf8),
    "archetype_id": pl.Utf8,
    "locus_score": pl.Float64,
    "regulator_domain_score": pl.Float64,
    "sensory_domain_score": pl.Float64,
    "proximity_score": pl.Float64,
    "archetype_conservation_score": pl.Float64,
    "enrichment_score": pl.Float64,
    "taxonomic_breadth_score": pl.Float64,
    "phylogenetic_profile_score": pl.Float64,
    "structural_plausibility_score": pl.Float64,
    "candidate_score": pl.Float64,
    "candidate_score_q": pl.Float64,
    "regulation_posterior": pl.Float64,
    "regulation_posterior_hdi_low": pl.Float64,
    "regulation_posterior_hdi_high": pl.Float64,
    "posterior_evidence_model": pl.Utf8,
    "rationale": pl.Utf8,
}


def _position(relative_index: int) -> str:
    if relative_index < 0:
        return "upstream"
    if relative_index > 0:
        return "downstream"
    return "internal"


def _distance_nt(row: dict[str, object]) -> int:
    return min(
        abs(int(cast(int, row["relative_start"]))),
        abs(int(cast(int, row["relative_stop"]))),
    )


def _dna_binding_domains(pfam_ids: list[str]) -> list[str]:
    return [pfam_id for pfam_id in pfam_ids if pfam_id in DNA_BINDING_PFAMS]


def _proximity_score(relative_index: int) -> float:
    return 1.0 / (1.0 + abs(relative_index))


def _enrichment_lookup(
    enrichment: pl.DataFrame | None,
) -> dict[tuple[str, str], tuple[float, float]]:
    if enrichment is None or enrichment.is_empty():
        return {}
    lookup: dict[tuple[str, str], tuple[float, float]] = {}
    for row in enrichment.iter_rows(named=True):
        odds_ratio = float(row["odds_ratio"])
        q_value = float(row["q_value"])
        enrichment_score = 0.0
        if q_value <= 0.05 and odds_ratio > 1.0:
            enrichment_score = math.log2(odds_ratio)
        lookup[(str(row["feature_type"]), str(row["feature_name"]))] = (
            enrichment_score,
            q_value,
        )
    return lookup


def _candidate_enrichment(
    row: dict[str, object],
    lookup: dict[tuple[str, str], tuple[float, float]],
) -> tuple[float, float | None]:
    scores: list[tuple[float, float]] = []
    regulator_match = lookup.get(("regulator_class", str(row["regulator_class"])))
    if regulator_match is not None:
        scores.append(regulator_match)
    for domain in cast(list[Any], row["sensory_domains"]):
        domain_match = lookup.get(("sensory_domain", str(domain)))
        if domain_match is not None:
            scores.append(domain_match)
    if not scores:
        return 0.0, None
    best_score = max(score for score, _ in scores)
    best_q = min(q_value for _, q_value in scores)
    return best_score, best_q


def _score_total(row: dict[str, object], scoring: ScoringConfig) -> float:
    weights = scoring.candidate_score_weights
    structural_score = row["structural_plausibility_score"]
    structural_component = (
        0.0 if structural_score is None else float(cast(float, structural_score))
    )
    return (
        float(cast(float, row["locus_score"])) * weights["locus"]
        + float(cast(float, row["regulator_domain_score"]))
        * weights["regulator_domain"]
        + float(cast(float, row["sensory_domain_score"])) * weights["sensory_domain"]
        + float(cast(float, row["proximity_score"])) * weights["proximity"]
        + float(cast(float, row["archetype_conservation_score"]))
        * weights["archetype_conservation"]
        + float(cast(float, row["enrichment_score"])) * weights["enrichment"]
        + float(cast(float, row["taxonomic_breadth_score"]))
        * weights.get("taxonomic_breadth", 0.0)
        + float(cast(float, row.get("phylogenetic_profile_score", 0.0)))
        * weights.get("phylogenetic_profile", 0.0)
        + structural_component * weights["structural_plausibility"]
    )


def score_candidates(
    loci: pl.DataFrame,
    genes: pl.DataFrame,
    scoring: ScoringConfig,
    enrichment: pl.DataFrame | None = None,
) -> pl.DataFrame:
    """Extract and score regulator candidates near scored loci."""

    locus_lookup = {str(row["locus_id"]): row for row in loci.iter_rows(named=True)}
    enrichment_scores = _enrichment_lookup(enrichment)
    rows: list[dict[str, object]] = []

    for gene in genes.filter(pl.col("is_regulator_candidate")).iter_rows(named=True):
        locus = locus_lookup[str(gene["locus_id"])]
        pfam_ids = [str(pfam_id) for pfam_id in cast(list[Any], gene["pfam_ids"])]
        sensory_domains = [
            str(domain) for domain in cast(list[Any], gene["sensory_domains"])
        ]
        enrichment_score, q_value = _candidate_enrichment(gene, enrichment_scores)
        row: dict[str, object] = {
            "candidate_id": f"{gene['locus_id']}::{gene['gene_accession']}",
            "analyte": locus["analyte"],
            "locus_id": gene["locus_id"],
            "gene_accession": gene["gene_accession"],
            "organism": locus["organism"],
            "cluster_id": locus["cluster_id"],
            "relative_index": gene["relative_index"],
            "distance_nt": _distance_nt(gene),
            "position": _position(int(gene["relative_index"])),
            "strand": gene["strand"],
            "regulator_class": gene["regulator_class"],
            "dna_binding_domains": _dna_binding_domains(pfam_ids),
            "sensory_domains": sensory_domains,
            "pfam_ids": pfam_ids,
            "interpro_ids": [
                str(item) for item in cast(list[Any], gene["interpro_ids"])
            ],
            "archetype_id": None,
            "locus_score": locus["locus_score"],
            "regulator_domain_score": 1.0,
            "sensory_domain_score": min(float(len(sensory_domains)), 2.0),
            "proximity_score": _proximity_score(int(gene["relative_index"])),
            "archetype_conservation_score": 0.0,
            "enrichment_score": enrichment_score,
            "taxonomic_breadth_score": 0.0,
            "phylogenetic_profile_score": 0.0,
            "structural_plausibility_score": None,
            "candidate_score": 0.0,
            "candidate_score_q": q_value,
            "regulation_posterior": None,
            "regulation_posterior_hdi_low": None,
            "regulation_posterior_hdi_high": None,
            "posterior_evidence_model": None,
            "rationale": "",
        }
        row["candidate_score"] = _score_total(row, scoring)
        row["rationale"] = (
            f"{row['regulator_class']} regulator {row['position']} of "
            f"{locus['anchor_family']} locus; "
            f"{len(sensory_domains)} sensory domain(s); "
            f"candidate score {float(cast(float, row['candidate_score'])):.3f}"
        )
        rows.append(row)

    if not rows:
        return validate(
            pl.DataFrame(schema=CANDIDATE_SCHEMA),
            RegulatorCandidatesSchema,
        )

    candidates = pl.DataFrame(rows, schema_overrides=CANDIDATE_SCHEMA).sort(
        ["candidate_score", "candidate_id"],
        descending=[True, False],
    )
    return validate(candidates, RegulatorCandidatesSchema)
