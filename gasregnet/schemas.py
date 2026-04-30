"""Pandera schemas for cross-module tables."""

from __future__ import annotations

from typing import Any, cast

import pandera.polars as pa
import polars as pl

from gasregnet.errors import SchemaError

DataFrameSchema = Any

ANALYTES = ["CO", "CN"]
LOCUS_CONFIDENCE = ["high", "medium", "low", "control"]
FUNCTIONAL_CLASSES = ["anchor", "regulator", "metabolic", "transporter", "unknown"]
REGULATOR_CLASSES = [
    "one_component",
    "two_component_rr",
    "two_component_hk",
    "sigma",
    "anti_sigma",
    "antiterminator",
    "none",
]
ANCHOR_EVIDENCE_TYPES = ["term_scan", "diamond", "hmmer"]
SENSOR_ROLES = ["anchor", "regulator", "sensor", "both", "accessory", "none"]


def _column(dtype: Any, **kwargs: Any) -> pa.Column:
    return pa.Column(dtype, **kwargs)


def _schema(columns: dict[str, pa.Column]) -> pa.DataFrameSchema:
    return pa.DataFrameSchema(columns, strict=True, coerce=False)


AnchorHitsSchema = _schema(
    {
        "dataset_name": _column(pl.Utf8),
        "analyte": _column(pl.Utf8, checks=pa.Check.isin(ANALYTES)),
        "anchor_family": _column(pl.Utf8),
        "protein_accession": _column(pl.Utf8),
        "locus_tag": _column(pl.Utf8),
        "gene": _column(pl.Utf8),
        "product": _column(pl.Utf8),
        "bitscore": _column(pl.Float64, nullable=True),
        "e_value": _column(pl.Float64, nullable=True),
        "identity": _column(pl.Float64, nullable=True),
        "coverage": _column(pl.Float64, nullable=True),
        "evidence_type": _column(
            pl.Utf8,
            checks=pa.Check.isin(ANCHOR_EVIDENCE_TYPES),
        ),
    },
)


LociSchema = _schema(
    {
        "locus_id": _column(pl.Utf8),
        "analyte": _column(pl.Utf8, checks=pa.Check.isin(ANALYTES)),
        "anchor_accession": _column(pl.Utf8),
        "anchor_family": _column(pl.Utf8),
        "organism": _column(pl.Utf8),
        "taxon_id": _column(pl.Int64),
        "cluster_id": _column(pl.Int32),
        "contig_id": _column(pl.Utf8),
        "window_size": _column(pl.Int32),
        "is_boundary_truncated": _column(pl.Boolean),
        "marker_genes_present": _column(pl.List(pl.Utf8)),
        "accessory_genes_present": _column(pl.List(pl.Utf8)),
        "locus_score": _column(pl.Float64),
        "locus_confidence": _column(
            pl.Utf8,
            checks=pa.Check.isin(LOCUS_CONFIDENCE),
        ),
        "taxonomic_context_score": _column(pl.Float64),
        "operon_integrity_score": _column(pl.Float64),
        "created_at": _column(pl.Datetime("us")),
        "phylum": _column(pl.Utf8, required=False),
        "class": _column(pl.Utf8, required=False),
        "order": _column(pl.Utf8, required=False),
        "family": _column(pl.Utf8, required=False),
        "genus": _column(pl.Utf8, required=False),
        "conservation_across_taxa_score": _column(pl.Float64, required=False),
    },
)

GenesSchema = _schema(
    {
        "locus_id": _column(pl.Utf8),
        "gene_accession": _column(pl.Utf8),
        "relative_index": _column(pl.Int32),
        "relative_start": _column(pl.Int64),
        "relative_stop": _column(pl.Int64),
        "strand": _column(pl.Utf8, checks=pa.Check.isin(["+", "-"])),
        "product_description": _column(pl.Utf8),
        "pfam_ids": _column(pl.List(pl.Utf8)),
        "pfam_descriptions": _column(pl.List(pl.Utf8)),
        "interpro_ids": _column(pl.List(pl.Utf8)),
        "interpro_descriptions": _column(pl.List(pl.Utf8)),
        "functional_class": _column(
            pl.Utf8,
            checks=pa.Check.isin(FUNCTIONAL_CLASSES),
        ),
        "regulator_class": _column(
            pl.Utf8,
            checks=pa.Check.isin(REGULATOR_CLASSES),
        ),
        "sensory_domains": _column(pl.List(pl.Utf8)),
        "is_anchor": _column(pl.Boolean),
        "is_regulator_candidate": _column(pl.Boolean),
        "sensor_role": _column(
            pl.Utf8,
            checks=pa.Check.isin(SENSOR_ROLES),
            required=False,
        ),
    },
)

SensorRegulatorPairsSchema = _schema(
    {
        "pair_id": _column(pl.Utf8),
        "analyte": _column(pl.Utf8, checks=pa.Check.isin(ANALYTES)),
        "locus_id": _column(pl.Utf8),
        "hk_gene_accession": _column(pl.Utf8),
        "rr_gene_accession": _column(pl.Utf8),
        "intergenic_distance_nt": _column(pl.Int64),
        "co_strand": _column(pl.Boolean),
        "hk_sensor_domains": _column(pl.List(pl.Utf8)),
        "hk_sensory_chemistries": _column(pl.List(pl.Utf8)),
        "rr_dna_binding_domains": _column(pl.List(pl.Utf8)),
        "pair_score": _column(pl.Float64),
    },
)

RegulatorCandidatesSchema = _schema(
    {
        "candidate_id": _column(pl.Utf8),
        "analyte": _column(pl.Utf8, checks=pa.Check.isin(ANALYTES)),
        "locus_id": _column(pl.Utf8),
        "gene_accession": _column(pl.Utf8),
        "organism": _column(pl.Utf8),
        "cluster_id": _column(pl.Int32),
        "relative_index": _column(pl.Int32),
        "distance_nt": _column(pl.Int64),
        "position": _column(
            pl.Utf8,
            checks=pa.Check.isin(["upstream", "downstream", "internal"]),
        ),
        "strand": _column(pl.Utf8, checks=pa.Check.isin(["+", "-"])),
        "regulator_class": _column(
            pl.Utf8,
            checks=pa.Check.isin(REGULATOR_CLASSES),
        ),
        "dna_binding_domains": _column(pl.List(pl.Utf8)),
        "sensory_domains": _column(pl.List(pl.Utf8)),
        "pfam_ids": _column(pl.List(pl.Utf8)),
        "interpro_ids": _column(pl.List(pl.Utf8)),
        "archetype_id": _column(pl.Utf8, nullable=True),
        "locus_score": _column(pl.Float64),
        "regulator_domain_score": _column(pl.Float64),
        "sensory_domain_score": _column(pl.Float64),
        "proximity_score": _column(pl.Float64),
        "archetype_conservation_score": _column(pl.Float64),
        "enrichment_score": _column(pl.Float64),
        "taxonomic_breadth_score": _column(pl.Float64),
        "structural_plausibility_score": _column(pl.Float64, nullable=True),
        "candidate_score": _column(pl.Float64),
        "candidate_score_q": _column(
            pl.Float64,
            checks=pa.Check.in_range(0.0, 1.0),
            nullable=True,
        ),
        "rationale": _column(pl.Utf8),
    },
)

ArchetypesSchema = _schema(
    {
        "archetype_id": _column(pl.Utf8),
        "analyte": _column(pl.Utf8, checks=pa.Check.isin(ANALYTES)),
        "cluster_id": _column(pl.Int32),
        "architecture_string": _column(pl.Utf8),
        "n_loci": _column(pl.Int32),
        "n_taxa": _column(pl.Int32),
        "representative_locus_id": _column(pl.Utf8),
        "dominant_regulator_class": _column(pl.Utf8),
        "dominant_anchor_structure": _column(pl.Utf8),
        "mean_locus_score": _column(pl.Float64),
        "mean_candidate_score": _column(pl.Float64),
    },
)

EnrichmentResultsSchema = _schema(
    {
        "analyte": _column(pl.Utf8, checks=pa.Check.isin(ANALYTES)),
        "feature_type": _column(
            pl.Utf8,
            checks=pa.Check.isin(
                [
                    "regulator_family",
                    "sensory_domain",
                    "regulator_class",
                    "archetype",
                ],
            ),
        ),
        "feature_name": _column(pl.Utf8),
        "case_definition": _column(pl.Utf8),
        "control_definition": _column(pl.Utf8),
        "n_case_with_feature": _column(pl.Int64),
        "n_case_without_feature": _column(pl.Int64),
        "n_control_with_feature": _column(pl.Int64),
        "n_control_without_feature": _column(pl.Int64),
        "odds_ratio": _column(pl.Float64),
        "p_value": _column(pl.Float64, checks=pa.Check.in_range(0.0, 1.0)),
        "q_value": _column(pl.Float64, checks=pa.Check.in_range(0.0, 1.0)),
        "interpretation": _column(pl.Utf8),
    },
)

BenchmarkSchema = _schema(
    {
        "benchmark_id": _column(pl.Utf8),
        "analyte": _column(
            pl.Utf8,
            checks=pa.Check.isin(["CO", "CN", "negative_control"]),
        ),
        "protein_name": _column(pl.Utf8),
        "uniprot_accession": _column(pl.Utf8),
        "organism": _column(pl.Utf8),
        "taxon_id": _column(pl.Int64),
        "anchor_family": _column(pl.Utf8, nullable=True),
        "expected_regulator_class": _column(pl.Utf8),
        "expected_sensory_domains": _column(pl.List(pl.Utf8)),
        "sensing_evidence_class": _column(
            pl.Utf8,
            checks=pa.Check.isin(["direct", "indirect", "inferred", "none"]),
        ),
        "pmid": _column(pl.List(pl.Utf8)),
        "notes": _column(pl.Utf8),
        "first_publication": _column(pl.Utf8, required=False),
        "verify_pmid": _column(pl.Boolean, required=False),
    },
)


def validate(df: pl.DataFrame, schema: DataFrameSchema) -> pl.DataFrame:
    """Validate a Polars DataFrame and raise a project-specific schema error."""

    try:
        validated = cast(pl.DataFrame, schema.validate(df))
    except Exception as exc:
        raise SchemaError(str(exc)) from exc

    if schema is GenesSchema:
        invalid_anchor_rows = validated.filter(
            (pl.col("relative_index") == 0) & (pl.col("is_anchor").not_()),
        )
        if invalid_anchor_rows.height > 0:
            raise SchemaError("relative_index == 0 requires is_anchor == True")
    return validated
