"""V2 gene role assignment for anchors, regulators, and sensors."""

from __future__ import annotations

from typing import Any

import polars as pl

from gasregnet.annotation.regulators import _classify_regulator
from gasregnet.config import (
    RegulatorFamilyEntry,
    SensoryDomainEntry,
    SensoryPairedEvidenceRule,
)
from gasregnet.schemas import GenesSchema, SensorRegulatorPairsSchema, validate

GENES_WITH_ROLES_SCHEMA_OVERRIDES: dict[str, Any] = {
    "relative_index": pl.Int32,
    "relative_start": pl.Int64,
    "relative_stop": pl.Int64,
    "pfam_ids": pl.List(pl.Utf8),
    "pfam_descriptions": pl.List(pl.Utf8),
    "interpro_ids": pl.List(pl.Utf8),
    "interpro_descriptions": pl.List(pl.Utf8),
    "sensory_domains": pl.List(pl.Utf8),
    "is_anchor": pl.Boolean,
    "is_regulator_candidate": pl.Boolean,
    "sensor_role": pl.Utf8,
}
PAIR_SCHEMA: dict[str, Any] = {
    "pair_id": pl.Utf8,
    "analyte": pl.Utf8,
    "locus_id": pl.Utf8,
    "hk_gene_accession": pl.Utf8,
    "rr_gene_accession": pl.Utf8,
    "intergenic_distance_nt": pl.Int64,
    "co_strand": pl.Boolean,
    "hk_sensor_domains": pl.List(pl.Utf8),
    "hk_sensory_chemistries": pl.List(pl.Utf8),
    "rr_dna_binding_domains": pl.List(pl.Utf8),
    "pair_score": pl.Float64,
}
DNA_BINDING_PFAMS = {
    "PF00126",
    "PF00158",
    "PF00196",
    "PF00325",
    "PF00356",
    "PF00376",
    "PF00392",
    "PF00440",
    "PF01022",
    "PF01047",
    "PF01381",
    "PF02954",
    "PF04542",
    "PF04545",
    "PF12833",
}


def _sensory_catalog(
    entries: list[SensoryDomainEntry],
) -> dict[str, SensoryDomainEntry]:
    return {entry.pfam_id: entry for entry in entries}


def _sensor_chemistries(
    pfam_ids: set[str],
    catalog: dict[str, SensoryDomainEntry],
    paired_evidence_rules: list[SensoryPairedEvidenceRule],
) -> list[str]:
    chemistries: list[str] = [
        entry.chemistry
        for pfam_id in pfam_ids
        if (entry := catalog.get(pfam_id)) is not None
        and entry.role == "sensor"
        and entry.chemistry != "none"
    ]
    for rule in paired_evidence_rules:
        if set(rule.if_pfam_all).issubset(pfam_ids) and (
            not rule.if_co_pfam_any or bool(set(rule.if_co_pfam_any) & pfam_ids)
        ):
            chemistries.append(rule.rescore.chemistry)
    return sorted(set(chemistries))


def _has_sensor_evidence(
    pfam_ids: set[str],
    catalog: dict[str, SensoryDomainEntry],
    paired_evidence_rules: list[SensoryPairedEvidenceRule],
) -> bool:
    if _sensor_chemistries(pfam_ids, catalog, paired_evidence_rules):
        return True
    return any(
        (entry := catalog.get(pfam_id)) is not None and entry.role == "sensor"
        for pfam_id in pfam_ids
    )


def _has_annotation(row: dict[str, object]) -> bool:
    for column in ("pfam_ids", "interpro_ids", "product_description"):
        value = row[column]
        if isinstance(value, list) and value:
            return True
        if isinstance(value, str) and value:
            return True
    return False


def _sensor_role(
    *,
    is_anchor: bool,
    regulator_class: str,
    sensor_evidence: bool,
    has_annotation: bool,
) -> str:
    if is_anchor:
        return "anchor"
    regulator_evidence = regulator_class != "none"
    if regulator_class == "two_component_hk" and sensor_evidence:
        return "sensor"
    if sensor_evidence and regulator_evidence:
        return "both"
    if regulator_evidence:
        return "regulator"
    if sensor_evidence:
        return "sensor"
    if has_annotation:
        return "accessory"
    return "none"


def assign_sensor_roles(
    genes: pl.DataFrame,
    *,
    regulator_families: list[RegulatorFamilyEntry],
    sensory_domain_catalog: list[SensoryDomainEntry],
    paired_evidence_rules: list[SensoryPairedEvidenceRule],
) -> pl.DataFrame:
    """Add sensor_role and refined regulator_class columns to genes."""

    catalog = _sensory_catalog(sensory_domain_catalog)
    rows: list[dict[str, object]] = []
    for row in genes.iter_rows(named=True):
        updated = dict(row)
        pfam_ids = {str(pfam_id) for pfam_id in row["pfam_ids"]}
        regulator_class = _classify_regulator(pfam_ids, regulator_families)
        sensor_evidence = _has_sensor_evidence(
            pfam_ids,
            catalog,
            paired_evidence_rules,
        )
        role = _sensor_role(
            is_anchor=bool(row["is_anchor"]),
            regulator_class=regulator_class,
            sensor_evidence=sensor_evidence,
            has_annotation=_has_annotation(row) or bool(row["is_regulator_candidate"]),
        )
        if regulator_class == "none" and bool(row["is_regulator_candidate"]):
            regulator_class = str(row["regulator_class"])
            if role == "sensor":
                role = "both"
            elif role in {"accessory", "none"}:
                role = "regulator"
        is_regulator = role in {"regulator", "both"}
        updated["sensor_role"] = role
        updated["regulator_class"] = regulator_class
        updated["is_regulator_candidate"] = is_regulator
        if role == "anchor":
            updated["functional_class"] = "anchor"
        elif is_regulator:
            updated["functional_class"] = "regulator"
        elif role == "sensor":
            updated["functional_class"] = "unknown"
        rows.append(updated)
    return validate(
        pl.DataFrame(rows, schema_overrides=GENES_WITH_ROLES_SCHEMA_OVERRIDES),
        GenesSchema,
    )


def _intergenic_distance(left: dict[str, object], right: dict[str, object]) -> int:
    left_stop = int(str(left["relative_stop"]))
    right_start = int(str(right["relative_start"]))
    right_stop = int(str(right["relative_stop"]))
    left_start = int(str(left["relative_start"]))
    if left_stop <= right_start:
        return right_start - left_stop
    if right_stop <= left_start:
        return left_start - right_stop
    return 0


def build_sensor_regulator_pairs(
    genes: pl.DataFrame,
    loci: pl.DataFrame,
    *,
    sensory_domain_catalog: list[SensoryDomainEntry],
    paired_evidence_rules: list[SensoryPairedEvidenceRule],
) -> pl.DataFrame:
    """Link same-locus histidine kinase sensors to response regulators."""

    catalog = _sensory_catalog(sensory_domain_catalog)
    analyte_by_locus = {
        str(row["locus_id"]): str(row["analyte"]) for row in loci.iter_rows(named=True)
    }
    rows: list[dict[str, object]] = []
    for partition in genes.partition_by("locus_id", maintain_order=True):
        locus_id = str(partition["locus_id"].item(0))
        locus_genes = list(partition.iter_rows(named=True))
        hks = [
            row
            for row in locus_genes
            if row["regulator_class"] == "two_component_hk"
            and row.get("sensor_role") in {"sensor", "both"}
        ]
        rrs = [
            row
            for row in locus_genes
            if row["regulator_class"] == "two_component_rr"
            and row.get("sensor_role") in {"regulator", "both"}
        ]
        for hk in hks:
            hk_pfams = {str(pfam_id) for pfam_id in hk["pfam_ids"]}
            hk_chemistries = _sensor_chemistries(
                hk_pfams,
                catalog,
                paired_evidence_rules,
            )
            hk_domains = [
                catalog[pfam_id].domain for pfam_id in hk_pfams if pfam_id in catalog
            ]
            for rr in rrs:
                rr_pfams = [str(pfam_id) for pfam_id in rr["pfam_ids"]]
                distance = _intergenic_distance(hk, rr)
                rows.append(
                    {
                        "pair_id": (
                            f"{locus_id}::{hk['gene_accession']}::"
                            f"{rr['gene_accession']}"
                        ),
                        "analyte": analyte_by_locus[str(locus_id)],
                        "locus_id": str(locus_id),
                        "hk_gene_accession": str(hk["gene_accession"]),
                        "rr_gene_accession": str(rr["gene_accession"]),
                        "intergenic_distance_nt": distance,
                        "co_strand": hk["strand"] == rr["strand"],
                        "hk_sensor_domains": sorted(set(hk_domains)),
                        "hk_sensory_chemistries": hk_chemistries,
                        "rr_dna_binding_domains": [
                            pfam_id
                            for pfam_id in rr_pfams
                            if pfam_id in DNA_BINDING_PFAMS
                        ],
                        "pair_score": 1.0 / (1.0 + max(distance, 0)),
                    },
                )
    if not rows:
        return validate(pl.DataFrame(schema=PAIR_SCHEMA), SensorRegulatorPairsSchema)
    return validate(
        pl.DataFrame(rows, schema_overrides=PAIR_SCHEMA),
        SensorRegulatorPairsSchema,
    )
