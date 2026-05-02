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


def _pfam_id_set(pfam_ids: object) -> set[str]:
    if pfam_ids is None:
        return set()
    if isinstance(pfam_ids, pl.Series):
        values: list[object] = pfam_ids.to_list()
    elif isinstance(pfam_ids, list | tuple | set):
        values = list(pfam_ids)
    else:
        return set()
    return {str(pfam_id) for pfam_id in values}


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
    with_inferred = genes.with_columns(
        pl.col("pfam_ids")
        .map_elements(
            lambda pfam_ids: _classify_regulator(
                _pfam_id_set(pfam_ids),
                regulator_families,
            ),
            return_dtype=pl.Utf8,
        )
        .alias("__inferred_regulator_class"),
        pl.col("pfam_ids")
        .map_elements(
            lambda pfam_ids: _has_sensor_evidence(
                _pfam_id_set(pfam_ids),
                catalog,
                paired_evidence_rules,
            ),
            return_dtype=pl.Boolean,
        )
        .alias("__sensor_evidence"),
        (
            (pl.col("pfam_ids").list.len() > 0)
            | (pl.col("interpro_ids").list.len() > 0)
            | (pl.col("product_description").str.len_bytes() > 0)
            | pl.col("is_regulator_candidate")
        ).alias("__has_annotation"),
    )
    with_roles = with_inferred.with_columns(
        pl.struct(
            [
                "is_anchor",
                "__inferred_regulator_class",
                "__sensor_evidence",
                "__has_annotation",
            ],
        )
        .map_elements(
            lambda row: _sensor_role(
                is_anchor=bool(row["is_anchor"]),
                regulator_class=str(row["__inferred_regulator_class"]),
                sensor_evidence=bool(row["__sensor_evidence"]),
                has_annotation=bool(row["__has_annotation"]),
            ),
            return_dtype=pl.Utf8,
        )
        .alias("__inferred_sensor_role"),
    )
    rescued_candidate = (pl.col("__inferred_regulator_class") == "none") & pl.col(
        "is_regulator_candidate",
    )
    assigned = (
        with_roles.with_columns(
            pl.when(rescued_candidate)
            .then(pl.col("regulator_class"))
            .otherwise(pl.col("__inferred_regulator_class"))
            .alias("regulator_class"),
            pl.when(rescued_candidate & (pl.col("__inferred_sensor_role") == "sensor"))
            .then(pl.lit("both"))
            .when(
                rescued_candidate
                & pl.col("__inferred_sensor_role").is_in(["accessory", "none"]),
            )
            .then(pl.lit("regulator"))
            .otherwise(pl.col("__inferred_sensor_role"))
            .alias("sensor_role"),
        )
        .with_columns(
            pl.col("sensor_role")
            .is_in(["regulator", "both"])
            .alias("is_regulator_candidate"),
        )
        .with_columns(
            pl.when(pl.col("sensor_role") == "anchor")
            .then(pl.lit("anchor"))
            .when(pl.col("is_regulator_candidate"))
            .then(pl.lit("regulator"))
            .when(pl.col("sensor_role") == "sensor")
            .then(pl.lit("unknown"))
            .otherwise(pl.col("functional_class"))
            .alias("functional_class"),
        )
        .drop(
            [
                "__inferred_regulator_class",
                "__sensor_evidence",
                "__has_annotation",
                "__inferred_sensor_role",
            ],
        )
    )
    return validate(
        pl.DataFrame(assigned, schema_overrides=GENES_WITH_ROLES_SCHEMA_OVERRIDES),
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
    hks = genes.filter(
        (pl.col("regulator_class") == "two_component_hk")
        & pl.col("sensor_role").is_in(["sensor", "both"]),
    ).select(
        "locus_id",
        pl.col("gene_accession").alias("hk_gene_accession"),
        pl.col("relative_start").alias("hk_relative_start"),
        pl.col("relative_stop").alias("hk_relative_stop"),
        pl.col("strand").alias("hk_strand"),
        pl.col("pfam_ids").alias("hk_pfam_ids"),
    )
    rrs = genes.filter(
        (pl.col("regulator_class") == "two_component_rr")
        & pl.col("sensor_role").is_in(["regulator", "both"]),
    ).select(
        "locus_id",
        pl.col("gene_accession").alias("rr_gene_accession"),
        pl.col("relative_start").alias("rr_relative_start"),
        pl.col("relative_stop").alias("rr_relative_stop"),
        pl.col("strand").alias("rr_strand"),
        pl.col("pfam_ids").alias("rr_pfam_ids"),
    )

    if hks.is_empty() or rrs.is_empty():
        return validate(pl.DataFrame(schema=PAIR_SCHEMA), SensorRegulatorPairsSchema)
    pairs = hks.join(rrs, on="locus_id", how="inner")
    if pairs.is_empty():
        return validate(pl.DataFrame(schema=PAIR_SCHEMA), SensorRegulatorPairsSchema)

    distance = (
        pl.when(pl.col("hk_relative_stop") <= pl.col("rr_relative_start"))
        .then(pl.col("rr_relative_start") - pl.col("hk_relative_stop"))
        .when(pl.col("rr_relative_stop") <= pl.col("hk_relative_start"))
        .then(pl.col("hk_relative_start") - pl.col("rr_relative_stop"))
        .otherwise(pl.lit(0))
    )
    scored = (
        pairs.join(loci.select("locus_id", "analyte"), on="locus_id", how="left")
        .with_columns(
            distance.cast(pl.Int64).alias("intergenic_distance_nt"),
            (pl.col("hk_strand") == pl.col("rr_strand")).alias("co_strand"),
            pl.col("hk_pfam_ids")
            .map_elements(
                lambda pfam_ids: sorted(
                    {
                        catalog[pfam_id].domain
                        for pfam_id in _pfam_id_set(pfam_ids)
                        if pfam_id in catalog
                    },
                ),
                return_dtype=pl.List(pl.Utf8),
            )
            .alias("hk_sensor_domains"),
            pl.col("hk_pfam_ids")
            .map_elements(
                lambda pfam_ids: _sensor_chemistries(
                    _pfam_id_set(pfam_ids),
                    catalog,
                    paired_evidence_rules,
                ),
                return_dtype=pl.List(pl.Utf8),
            )
            .alias("hk_sensory_chemistries"),
            pl.col("rr_pfam_ids")
            .map_elements(
                lambda pfam_ids: [
                    str(pfam_id)
                    for pfam_id in pfam_ids
                    if str(pfam_id) in DNA_BINDING_PFAMS
                ],
                return_dtype=pl.List(pl.Utf8),
            )
            .alias("rr_dna_binding_domains"),
        )
        .with_columns(
            pl.concat_str(
                ["locus_id", "hk_gene_accession", "rr_gene_accession"],
                separator="::",
            ).alias("pair_id"),
            (1.0 / (1.0 + pl.col("intergenic_distance_nt"))).alias("pair_score"),
        )
        .select(list(PAIR_SCHEMA))
    )
    return validate(
        pl.DataFrame(scored, schema_overrides=PAIR_SCHEMA),
        SensorRegulatorPairsSchema,
    )
