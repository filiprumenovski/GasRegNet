from __future__ import annotations

from collections.abc import Callable
from datetime import datetime

import polars as pl
import pytest

from gasregnet import schemas
from gasregnet.errors import SchemaError

FrameFactory = Callable[[], pl.DataFrame]


def loci_frame() -> pl.DataFrame:
    return pl.DataFrame(
        {
            "locus_id": ["CO_1_anchor"],
            "analyte": ["CO"],
            "anchor_accession": ["anchor"],
            "anchor_family": ["coxL"],
            "organism": ["Rhodospirillum rubrum"],
            "taxon_id": [1085],
            "cluster_id": pl.Series([1], dtype=pl.Int32),
            "contig_id": ["contig"],
            "window_size": pl.Series([10], dtype=pl.Int32),
            "is_boundary_truncated": [False],
            "marker_genes_present": [["coxL"]],
            "accessory_genes_present": [["coxM"]],
            "locus_score": [6.0],
            "locus_confidence": ["high"],
            "taxonomic_context_score": [0.5],
            "operon_integrity_score": [1.0],
            "created_at": pl.Series([datetime(2026, 4, 29)], dtype=pl.Datetime("us")),
        },
    )


def anchor_hits_frame() -> pl.DataFrame:
    return pl.DataFrame(
        {
            "dataset_name": ["ecoli_k12_mg1655"],
            "analyte": ["CN"],
            "anchor_family": ["cydA"],
            "protein_accession": ["NP_415261.2"],
            "locus_tag": ["b0733"],
            "gene": ["cydA"],
            "product": ["cytochrome bd-I ubiquinol oxidase subunit I"],
            "bitscore": [None],
            "e_value": [None],
            "identity": [None],
            "coverage": [None],
            "evidence_type": ["term_scan"],
        },
        schema_overrides={
            "bitscore": pl.Float64,
            "e_value": pl.Float64,
            "identity": pl.Float64,
            "coverage": pl.Float64,
        },
    )


def genes_frame() -> pl.DataFrame:
    return pl.DataFrame(
        {
            "locus_id": ["CO_1_anchor"],
            "gene_accession": ["anchor"],
            "relative_index": pl.Series([0], dtype=pl.Int32),
            "relative_start": [1],
            "relative_stop": [100],
            "strand": ["+"],
            "product_description": ["carbon monoxide dehydrogenase"],
            "pfam_ids": [["PF02738"]],
            "pfam_descriptions": [["coxL"]],
            "interpro_ids": [["IPR000000"]],
            "interpro_descriptions": [["domain"]],
            "functional_class": ["anchor"],
            "regulator_class": ["none"],
            "sensory_domains": pl.Series([[]], dtype=pl.List(pl.Utf8)),
            "is_anchor": [True],
            "is_regulator_candidate": [False],
        },
    )


def candidates_frame() -> pl.DataFrame:
    return pl.DataFrame(
        {
            "candidate_id": ["cand1"],
            "analyte": ["CO"],
            "locus_id": ["CO_1_anchor"],
            "gene_accession": ["reg1"],
            "organism": ["Rhodospirillum rubrum"],
            "cluster_id": pl.Series([1], dtype=pl.Int32),
            "relative_index": pl.Series([-1], dtype=pl.Int32),
            "distance_nt": [150],
            "position": ["upstream"],
            "strand": ["+"],
            "regulator_class": ["one_component"],
            "dna_binding_domains": [["HTH"]],
            "sensory_domains": [["PAS"]],
            "pfam_ids": [["PF00989"]],
            "interpro_ids": [["IPR000000"]],
            "archetype_id": ["arch1"],
            "locus_score": [6.0],
            "regulator_domain_score": [1.5],
            "sensory_domain_score": [2.0],
            "proximity_score": [1.0],
            "archetype_conservation_score": [1.5],
            "enrichment_score": [2.0],
            "taxonomic_breadth_score": [0.5],
            "structural_plausibility_score": [None],
            "candidate_score": [14.5],
            "candidate_score_q": [0.01],
            "rationale": ["near cox locus"],
        },
        schema_overrides={
            "structural_plausibility_score": pl.Float64,
            "candidate_score_q": pl.Float64,
        },
    )


def archetypes_frame() -> pl.DataFrame:
    return pl.DataFrame(
        {
            "archetype_id": ["arch1"],
            "analyte": ["CO"],
            "cluster_id": pl.Series([1], dtype=pl.Int32),
            "architecture_string": ["[-1:TF:PAS][0:coxL]"],
            "n_loci": pl.Series([1], dtype=pl.Int32),
            "n_taxa": pl.Series([1], dtype=pl.Int32),
            "representative_locus_id": ["CO_1_anchor"],
            "dominant_regulator_class": ["one_component"],
            "dominant_anchor_structure": ["coxL"],
            "mean_locus_score": [6.0],
            "mean_candidate_score": [14.5],
        },
    )


def enrichment_frame() -> pl.DataFrame:
    return pl.DataFrame(
        {
            "analyte": ["CO"],
            "feature_type": ["sensory_domain"],
            "feature_name": ["PAS"],
            "case_definition": ["high confidence CO"],
            "control_definition": ["matched non-gas oxidase"],
            "n_case_with_feature": [10],
            "n_case_without_feature": [5],
            "n_control_with_feature": [3],
            "n_control_without_feature": [30],
            "odds_ratio": [20.0],
            "p_value": [0.001],
            "q_value": [0.01],
            "interpretation": ["enriched"],
        },
    )


def enrichment_robustness_frame() -> pl.DataFrame:
    return pl.DataFrame(
        {
            "analyte": ["CO"],
            "feature_type": ["sensory_domain"],
            "feature_name": ["PAS"],
            "deduplication_policy": ["one_per_family"],
            "stratum_column": ["genus"],
            "test": ["cmh"],
            "q_value": [0.01],
            "p_value": [0.001],
            "odds_ratio": [20.0],
        },
    )


def benchmark_frame() -> pl.DataFrame:
    return pl.DataFrame(
        {
            "benchmark_id": ["co_CooA_Rrubrum"],
            "analyte": ["CO"],
            "protein_name": ["CooA"],
            "uniprot_accession": ["P19919"],
            "organism": ["Rhodospirillum rubrum"],
            "taxon_id": [1085],
            "anchor_family": ["coxL"],
            "expected_regulator_class": ["one_component"],
            "expected_sensory_domains": [["heme"]],
            "sensing_evidence_class": ["direct"],
            "pmid": [["00000000"]],
            "notes": ["placeholder"],
        },
    )


@pytest.mark.parametrize(
    ("schema", "factory"),
    [
        (schemas.AnchorHitsSchema, anchor_hits_frame),
        (schemas.LociSchema, loci_frame),
        (schemas.GenesSchema, genes_frame),
        (schemas.RegulatorCandidatesSchema, candidates_frame),
        (schemas.ArchetypesSchema, archetypes_frame),
        (schemas.EnrichmentResultsSchema, enrichment_frame),
        (schemas.EnrichmentRobustnessSchema, enrichment_robustness_frame),
        (schemas.BenchmarkSchema, benchmark_frame),
    ],
)
def test_schema_accepts_valid_frame(
    schema: schemas.DataFrameSchema,
    factory: FrameFactory,
) -> None:
    assert schemas.validate(factory(), schema).height == 1


@pytest.mark.parametrize(
    ("schema", "factory"),
    [
        (schemas.AnchorHitsSchema, anchor_hits_frame),
        (schemas.LociSchema, loci_frame),
        (schemas.GenesSchema, genes_frame),
        (schemas.RegulatorCandidatesSchema, candidates_frame),
        (schemas.ArchetypesSchema, archetypes_frame),
        (schemas.EnrichmentResultsSchema, enrichment_frame),
        (schemas.EnrichmentRobustnessSchema, enrichment_robustness_frame),
        (schemas.BenchmarkSchema, benchmark_frame),
    ],
)
def test_schema_rejects_missing_column(
    schema: schemas.DataFrameSchema,
    factory: FrameFactory,
) -> None:
    frame = factory()

    with pytest.raises(SchemaError):
        schemas.validate(frame.drop(frame.columns[0]), schema)


@pytest.mark.parametrize(
    ("schema", "factory"),
    [
        (schemas.AnchorHitsSchema, anchor_hits_frame),
        (schemas.LociSchema, loci_frame),
        (schemas.GenesSchema, genes_frame),
        (schemas.RegulatorCandidatesSchema, candidates_frame),
        (schemas.ArchetypesSchema, archetypes_frame),
        (schemas.EnrichmentResultsSchema, enrichment_frame),
        (schemas.EnrichmentRobustnessSchema, enrichment_robustness_frame),
        (schemas.BenchmarkSchema, benchmark_frame),
    ],
)
def test_schema_rejects_wrong_dtype(
    schema: schemas.DataFrameSchema,
    factory: FrameFactory,
) -> None:
    frame = factory().with_columns(pl.lit(1).alias("analyte"))

    with pytest.raises(SchemaError):
        schemas.validate(frame, schema)


def test_genes_schema_enforces_anchor_cross_field_rule() -> None:
    frame = genes_frame().with_columns(pl.lit(False).alias("is_anchor"))

    with pytest.raises(SchemaError, match="relative_index == 0"):
        schemas.validate(frame, schemas.GenesSchema)
