"""Archetype conservation scoring."""

from __future__ import annotations

from typing import Any, cast

import polars as pl

from gasregnet.archetypes.cluster import architecture_string
from gasregnet.config import ScoringConfig
from gasregnet.schemas import RegulatorCandidatesSchema, validate
from gasregnet.scoring.candidates import candidate_score_from_components

CANDIDATE_SCHEMA_OVERRIDES: dict[str, Any] = {
    "cluster_id": pl.Int32,
    "relative_index": pl.Int32,
    "distance_nt": pl.Int64,
    "dna_binding_domains": pl.List(pl.Utf8),
    "sensory_domains": pl.List(pl.Utf8),
    "primary_sensory_chemistry": pl.Utf8,
    "pfam_ids": pl.List(pl.Utf8),
    "interpro_ids": pl.List(pl.Utf8),
    "archetype_id": pl.Utf8,
    "phylogenetic_profile_score": pl.Float64,
    "structural_plausibility_score": pl.Float64,
    "candidate_score_q": pl.Float64,
    "regulation_logit_score": pl.Float64,
    "score_band_low": pl.Float64,
    "score_band_high": pl.Float64,
    "score_band_model": pl.Utf8,
}


def _clip(value: float) -> float:
    return max(0.0, min(1.0, value))


def _taxonomy(row: dict[str, object], column: str) -> str:
    value = row.get(column)
    if isinstance(value, str) and value:
        return value
    organism = str(row.get("organism", ""))
    if column == "genus" and organism:
        return organism.split()[0]
    if column == "family":
        taxon_id = row.get("taxon_id", "")
        return str(taxon_id) if taxon_id is not None else ""
    if column == "phylum":
        return str(row.get("analyte", ""))
    return ""


def _locus_architecture_lookup(
    candidates: pl.DataFrame,
    loci: pl.DataFrame,
) -> dict[str, str]:
    candidates_by_locus: dict[str, list[dict[str, object]]] = {}
    for candidate in candidates.iter_rows(named=True):
        candidates_by_locus.setdefault(str(candidate["locus_id"]), []).append(candidate)
    lookup: dict[str, str] = {}
    for locus in loci.iter_rows(named=True):
        locus_id = str(locus["locus_id"])
        lookup[locus_id] = architecture_string(
            locus,
            candidates_by_locus.get(locus_id, []),
        )
    return lookup


def _archetype_lookup(archetypes: pl.DataFrame) -> dict[tuple[str, str], str]:
    lookup: dict[tuple[str, str], str] = {}
    for row in archetypes.iter_rows(named=True):
        architecture = str(row["architecture_string"])
        archetype_id = str(row["archetype_id"])
        lookup[(str(row["analyte"]), architecture)] = archetype_id
        lookup.setdefault(("", architecture), archetype_id)
    return lookup


def compute_conservation_scores(
    candidates: pl.DataFrame,
    archetypes: pl.DataFrame,
    loci: pl.DataFrame,
    *,
    min_loci_per_archetype: int = 3,
    min_taxa_per_archetype: int = 3,
    scoring: ScoringConfig | None = None,
) -> pl.DataFrame:
    """Compute archetype conservation and taxonomic breadth per candidate."""

    candidates = validate(candidates, RegulatorCandidatesSchema)
    if candidates.is_empty() or archetypes.is_empty() or loci.is_empty():
        return candidates

    architecture_by_locus = _locus_architecture_lookup(candidates, loci)
    archetype_by_architecture = _archetype_lookup(archetypes)
    locus_by_id = {str(row["locus_id"]): row for row in loci.iter_rows(named=True)}
    loci_by_architecture: dict[str, dict[str, dict[str, object]]] = {}
    candidates_by_arch_pos: dict[tuple[str, str, int, str], int] = {}
    for candidate in candidates.iter_rows(named=True):
        locus_id = str(candidate["locus_id"])
        architecture = architecture_by_locus.get(locus_id, "")
        if not architecture:
            continue
        locus = locus_by_id[locus_id]
        loci_by_architecture.setdefault(architecture, {})[locus_id] = locus
        key = (
            str(candidate["analyte"]),
            architecture,
            int(cast(int, candidate["relative_index"])),
            str(candidate["regulator_class"]),
        )
        candidates_by_arch_pos[key] = candidates_by_arch_pos.get(key, 0) + 1

    rows: list[dict[str, object]] = []
    for candidate in candidates.iter_rows(named=True):
        updated = dict(candidate)
        locus_id = str(candidate["locus_id"])
        architecture = architecture_by_locus.get(locus_id, "")
        arch_loci = list(loci_by_architecture.get(architecture, {}).values())
        archetype_id = archetype_by_architecture.get(
            (str(candidate["analyte"]), architecture),
            archetype_by_architecture.get(
                ("", architecture),
                str(candidate["archetype_id"] or ""),
            ),
        )
        updated["archetype_id"] = archetype_id
        genus_count = len(
            {genus for locus in arch_loci if (genus := _taxonomy(locus, "genus"))},
        )
        family_count = len(
            {family for locus in arch_loci if (family := _taxonomy(locus, "family"))},
        )
        taxon_count = max(genus_count, family_count)
        if (
            len(arch_loci) < min_loci_per_archetype
            or taxon_count < min_taxa_per_archetype
        ):
            updated["archetype_conservation_score"] = 0.0
            updated["taxonomic_breadth_score"] = 0.0
        else:
            n_loci = float(len(arch_loci))
            genus_breadth = _clip(
                len({_taxonomy(locus, "genus") for locus in arch_loci}) / n_loci,
            )
            family_breadth = _clip(
                len({_taxonomy(locus, "family") for locus in arch_loci}) / n_loci,
            )
            position_key = (
                str(candidate["analyte"]),
                architecture,
                int(cast(int, candidate["relative_index"])),
                str(candidate["regulator_class"]),
            )
            position_conservation = _clip(
                candidates_by_arch_pos.get(position_key, 0) / n_loci,
            )
            updated["archetype_conservation_score"] = (
                genus_breadth * family_breadth * position_conservation
            ) ** (1.0 / 3.0)
            updated["taxonomic_breadth_score"] = _clip(
                len({_taxonomy(locus, "phylum") for locus in arch_loci}) / 30.0,
            )
        if scoring is not None:
            updated["candidate_score"] = candidate_score_from_components(
                updated,
                scoring,
            )
        rows.append(updated)

    if not rows:
        return candidates
    return validate(
        pl.DataFrame(rows, schema_overrides=CANDIDATE_SCHEMA_OVERRIDES),
        RegulatorCandidatesSchema,
    )
