"""Candidate regulator scoring."""

from __future__ import annotations

import math
from typing import Any, cast

import polars as pl

from gasregnet.config import (
    AnalyteConfig,
    ScoringConfig,
    SensoryDomainEntry,
    SensoryPairedEvidenceRule,
)
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
    "primary_sensory_chemistry": pl.Utf8,
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
    "regulation_logit_score": pl.Float64,
    "score_band_low": pl.Float64,
    "score_band_high": pl.Float64,
    "score_band_model": pl.Utf8,
    "rationale": pl.Utf8,
}
CANDIDATE_SCORE_COMPONENTS = {
    "locus": "locus_score",
    "regulator_domain": "regulator_domain_score",
    "sensory_domain": "sensory_domain_score",
    "proximity": "proximity_score",
    "archetype_conservation": "archetype_conservation_score",
    "enrichment": "enrichment_score",
    "taxonomic_breadth": "taxonomic_breadth_score",
    "phylogenetic_profile": "phylogenetic_profile_score",
    "structural_plausibility": "structural_plausibility_score",
    "operator_motif": "operator_motif_score",
    "embedding_similarity": "embedding_similarity_score",
    "foldseek_similarity": "foldseek_similarity_score",
    "paired_sensor_evidence": "paired_sensor_evidence_score",
    "conservation_across_taxa": "conservation_across_taxa_score",
}


def expected_chemistry_by_analyte(
    analytes: list[AnalyteConfig],
) -> dict[str, set[str]]:
    """Return expected sensory chemistries keyed by analyte name."""

    return {
        analyte.analyte: {str(value) for value in analyte.expected_sensory_chemistry}
        for analyte in analytes
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
    scoring: ScoringConfig,
) -> dict[tuple[str, str], tuple[float, float]]:
    if enrichment is None or enrichment.is_empty():
        return {}
    lookup: dict[tuple[str, str], tuple[float, float]] = {}
    score_cap = scoring.enrichment.score_cap
    for row in enrichment.iter_rows(named=True):
        odds_ratio = float(row["odds_ratio"])
        q_value = float(row["q_value"])
        enrichment_score = 0.0
        if q_value <= 0.05 and odds_ratio > 1.0:
            enrichment_score = min(math.log2(odds_ratio), score_cap)
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


def candidate_score_from_components(
    row: dict[str, object],
    scoring: ScoringConfig,
) -> float:
    """Compute the decomposable candidate score from configured components."""

    total = 0.0
    for weight_name, weight in scoring.candidate_score_weights.items():
        component_name = CANDIDATE_SCORE_COMPONENTS.get(weight_name)
        if component_name is None:
            if weight == 0.0:
                continue
            raise ValueError(f"candidate score weight has no component: {weight_name}")
        component = row.get(component_name, 0.0)
        component_value = 0.0 if component is None else float(cast(float, component))
        total += component_value * float(weight)
    return total


def _regulator_domain_score(row: dict[str, object]) -> float:
    pfam_ids = {str(pfam_id) for pfam_id in cast(list[Any], row["pfam_ids"])}
    regulator_class = str(row["regulator_class"])
    if regulator_class == "none":
        return 0.0
    if _dna_binding_domains(sorted(pfam_ids)):
        return 1.0
    if regulator_class in {"two_component_rr", "two_component_hk", "sigma54_activator"}:
        return 0.75
    return 0.25


def _catalog_by_domain(
    entries: list[SensoryDomainEntry] | None,
) -> dict[str, SensoryDomainEntry]:
    if entries is None:
        return {}
    return {entry.domain: entry for entry in entries}


def _catalog_by_pfam(
    entries: list[SensoryDomainEntry] | None,
) -> dict[str, SensoryDomainEntry]:
    if entries is None:
        return {}
    return {entry.pfam_id: entry for entry in entries}


def _evidence_terms(*values: object) -> set[str]:
    terms: set[str] = set()
    for value in values:
        if value is None:
            continue
        if isinstance(value, str):
            terms.add(value.lower())
            continue
        if isinstance(value, list | tuple | set):
            terms.update(str(item).lower() for item in value)
    return terms


def _paired_evidence_rule_matches(
    rule: SensoryPairedEvidenceRule,
    *,
    pfam_ids: set[str],
    evidence_terms: set[str],
    organism_kingdom: str | None,
) -> bool:
    if not set(rule.if_pfam_all).issubset(pfam_ids):
        return False
    if rule.if_co_pfam_any and not bool(set(rule.if_co_pfam_any) & pfam_ids):
        return False
    if rule.if_motif_any and not bool(
        {motif.lower() for motif in rule.if_motif_any} & evidence_terms,
    ):
        return False
    if (
        rule.if_organism_kingdom is not None
        and (organism_kingdom or "").lower() != rule.if_organism_kingdom.lower()
    ):
        return False
    return True


def _sensor_chemistries(
    *,
    domains: list[str],
    pfam_ids: set[str],
    sensory_domain_catalog: list[SensoryDomainEntry] | None,
    paired_evidence_rules: list[SensoryPairedEvidenceRule] | None,
    evidence_terms: set[str] | None = None,
    organism_kingdom: str | None = None,
) -> list[str]:
    by_domain = _catalog_by_domain(sensory_domain_catalog)
    by_pfam = _catalog_by_pfam(sensory_domain_catalog)
    chemistries: list[str] = []
    for domain in domains:
        entry = by_domain.get(domain)
        if entry is not None and entry.role == "sensor" and entry.chemistry != "none":
            chemistries.append(entry.chemistry)
    for pfam_id in pfam_ids:
        entry = by_pfam.get(pfam_id)
        if entry is not None and entry.role == "sensor" and entry.chemistry != "none":
            chemistries.append(entry.chemistry)
    for rule in paired_evidence_rules or []:
        if _paired_evidence_rule_matches(
            rule,
            pfam_ids=pfam_ids,
            evidence_terms=evidence_terms or set(),
            organism_kingdom=organism_kingdom,
        ):
            chemistries.append(rule.rescore.chemistry)
    return sorted(set(chemistries))


def _sensory_domain_score(
    *,
    analyte: str,
    sensory_domains: list[str],
    pfam_ids: set[str],
    sensory_domain_catalog: list[SensoryDomainEntry] | None,
    paired_evidence_rules: list[SensoryPairedEvidenceRule] | None,
    expected_chemistry_by_analyte: dict[str, set[str]] | None,
    evidence_terms: set[str] | None = None,
    organism_kingdom: str | None = None,
) -> tuple[float, str]:
    if sensory_domain_catalog is None or expected_chemistry_by_analyte is None:
        primary = sensory_domains[0] if sensory_domains else "none"
        return min(float(len(sensory_domains)), 2.0), primary
    chemistries = _sensor_chemistries(
        domains=sensory_domains,
        pfam_ids=pfam_ids,
        sensory_domain_catalog=sensory_domain_catalog,
        paired_evidence_rules=paired_evidence_rules,
        evidence_terms=evidence_terms,
        organism_kingdom=organism_kingdom,
    )
    if not chemistries:
        return 0.0, "none"
    expected = expected_chemistry_by_analyte.get(analyte, set())
    matches = [chemistry for chemistry in chemistries if chemistry in expected]
    primary = matches[0] if matches else chemistries[0]
    score = float(len(matches)) if matches else 0.25
    return min(score, 2.0), primary


def score_candidates(
    loci: pl.DataFrame,
    genes: pl.DataFrame,
    scoring: ScoringConfig,
    enrichment: pl.DataFrame | None = None,
    *,
    sensory_domain_catalog: list[SensoryDomainEntry] | None = None,
    paired_evidence_rules: list[SensoryPairedEvidenceRule] | None = None,
    expected_chemistry_by_analyte: dict[str, set[str]] | None = None,
) -> pl.DataFrame:
    """Extract and score regulator candidates near scored loci."""

    locus_lookup = {str(row["locus_id"]): row for row in loci.iter_rows(named=True)}
    enrichment_scores = _enrichment_lookup(enrichment, scoring)
    rows: list[dict[str, object]] = []

    for gene in genes.filter(pl.col("is_regulator_candidate")).iter_rows(named=True):
        locus = locus_lookup[str(gene["locus_id"])]
        pfam_ids = [str(pfam_id) for pfam_id in cast(list[Any], gene["pfam_ids"])]
        sensory_domains = [
            str(domain) for domain in cast(list[Any], gene["sensory_domains"])
        ]
        sensory_score, primary_chemistry = _sensory_domain_score(
            analyte=str(locus["analyte"]),
            sensory_domains=sensory_domains,
            pfam_ids=set(pfam_ids),
            sensory_domain_catalog=sensory_domain_catalog,
            paired_evidence_rules=paired_evidence_rules,
            expected_chemistry_by_analyte=expected_chemistry_by_analyte,
            evidence_terms=_evidence_terms(
                sensory_domains,
                gene.get("pfam_descriptions"),
                gene.get("interpro_descriptions"),
                gene.get("product_description"),
            ),
            organism_kingdom=str(locus.get("superkingdom") or ""),
        )
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
            "primary_sensory_chemistry": primary_chemistry,
            "pfam_ids": pfam_ids,
            "interpro_ids": [
                str(item) for item in cast(list[Any], gene["interpro_ids"])
            ],
            "archetype_id": None,
            "locus_score": locus["locus_score"],
            "regulator_domain_score": _regulator_domain_score(gene),
            "sensory_domain_score": sensory_score,
            "proximity_score": _proximity_score(int(gene["relative_index"])),
            "archetype_conservation_score": 0.0,
            "enrichment_score": enrichment_score,
            "taxonomic_breadth_score": 0.0,
            "phylogenetic_profile_score": 0.0,
            "structural_plausibility_score": None,
            "candidate_score": 0.0,
            "candidate_score_q": q_value,
            "regulation_logit_score": None,
            "score_band_low": None,
            "score_band_high": None,
            "score_band_model": None,
            "rationale": "",
        }
        row["candidate_score"] = candidate_score_from_components(row, scoring)
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
