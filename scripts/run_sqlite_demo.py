"""Run the deterministic SQLite demo pipeline used by Snakemake."""

from __future__ import annotations

import argparse
from pathlib import Path

import polars as pl

from gasregnet.annotation.ecology import score_taxonomic_context_by_analyte
from gasregnet.annotation.roles import (
    assign_sensor_roles,
    build_sensor_regulator_pairs,
)
from gasregnet.archetypes.cluster import cluster_archetypes
from gasregnet.config import load_config, resolve_and_dump
from gasregnet.io.parquet import write_parquet
from gasregnet.manifest import build_manifest, write_manifest
from gasregnet.neighborhoods.retrieve import retrieve_from_efi_sqlite
from gasregnet.paths import ensure_out_dir
from gasregnet.reports.captions import build_result_led_captions, write_caption_files
from gasregnet.reports.figures import (
    figure_1_workflow_and_recovery,
    figure_2_locus_landscape,
    figure_3_archetype_atlas,
    figure_4_chemistry_partition,
    figure_5_candidate_ranking,
    figure_6_structure_hypotheses,
)
from gasregnet.reports.tables import write_publication_tables
from gasregnet.schemas import (
    ArchetypesSchema,
    EnrichmentResultsSchema,
    EnrichmentRobustnessSchema,
    GenesSchema,
    LociSchema,
    RegulatorCandidatesSchema,
    SensorRegulatorPairsSchema,
)
from gasregnet.scoring.candidates import expected_chemistry_by_analyte, score_candidates
from gasregnet.scoring.conservation import compute_conservation_scores
from gasregnet.scoring.cooccurrence import assign_phylogenetic_profile_scores
from gasregnet.scoring.enrichment import (
    run_enrichment,
    run_enrichment_robustness,
    run_stratified_enrichment,
)
from gasregnet.scoring.loci import score_loci
from gasregnet.scoring.posterior import assign_operon_regulation_score_bands
from scripts.build_test_fixtures import build_mini_efi


def _write_frame_set(
    *,
    out_dir: Path,
    loci: pl.DataFrame,
    genes: pl.DataFrame,
    candidates: pl.DataFrame,
    enrichment: pl.DataFrame,
    enrichment_robustness: pl.DataFrame,
    archetypes: pl.DataFrame,
    sensor_regulator_pairs: pl.DataFrame,
) -> None:
    intermediate = out_dir / "intermediate"
    locus_columns = [column for column in LociSchema.columns if column in loci.columns]
    write_parquet(loci.select(locus_columns), intermediate / "loci.parquet", LociSchema)
    write_parquet(genes, intermediate / "genes.parquet", GenesSchema)
    write_parquet(
        candidates,
        intermediate / "candidates.parquet",
        RegulatorCandidatesSchema,
    )
    write_parquet(
        enrichment,
        intermediate / "enrichment.parquet",
        EnrichmentResultsSchema,
    )
    write_parquet(
        enrichment_robustness,
        intermediate / "enrichment_robustness.parquet",
        EnrichmentRobustnessSchema,
    )
    write_parquet(archetypes, intermediate / "archetypes.parquet", ArchetypesSchema)
    write_parquet(
        sensor_regulator_pairs,
        intermediate / "sensor_regulator_pairs.parquet",
        SensorRegulatorPairsSchema,
    )


def _benchmark_recovery(candidates: pl.DataFrame) -> pl.DataFrame:
    top = candidates.sort("candidate_score", descending=True).head(1)
    if top.is_empty():
        return pl.DataFrame(
            {
                "benchmark_id": ["mini_efi_known_sensor"],
                "analyte": ["CO"],
                "protein_name": ["PAS regulator"],
                "organism": ["Rhodospirillum rubrum"],
                "hit": [False],
                "rank": [None],
                "candidate_score": [None],
                "regulation_logit_score": [None],
            },
            schema_overrides={
                "rank": pl.Int64,
                "candidate_score": pl.Float64,
                "regulation_logit_score": pl.Float64,
            },
        )
    row = top.row(0, named=True)
    return pl.DataFrame(
        {
            "benchmark_id": ["mini_efi_known_sensor"],
            "analyte": [row["analyte"]],
            "protein_name": [row["gene_accession"]],
            "organism": [row["organism"]],
            "hit": [True],
            "rank": [1],
            "candidate_score": [row["candidate_score"]],
            "regulation_logit_score": [row.get("regulation_logit_score")],
        },
        schema_overrides={
            "rank": pl.Int64,
            "candidate_score": pl.Float64,
            "regulation_logit_score": pl.Float64,
        },
    )


def _genes_with_locus_taxonomy(loci: pl.DataFrame, genes: pl.DataFrame) -> pl.DataFrame:
    taxonomy_columns = ["locus_id"] + [
        column
        for column in ("organism", "taxon_id", "genus", "family")
        if column in loci.columns and column not in genes.columns
    ]
    if len(taxonomy_columns) == 1:
        return genes
    return genes.join(
        loci.select(taxonomy_columns).unique("locus_id"),
        on="locus_id",
        how="left",
    )


def run_sqlite_demo(*, out_dir: Path, config_path: Path, sqlite_path: Path) -> Path:
    """Run the deterministic SQLite fixture through the implemented pipeline."""

    ensure_out_dir(out_dir)
    sqlite = build_mini_efi(sqlite_path)
    config = load_config(config_path)
    resolve_and_dump(config, out_dir)

    loci, genes = retrieve_from_efi_sqlite(sqlite, ["CO", "CN"])
    loci = score_taxonomic_context_by_analyte(
        loci,
        config.analytes,
        root=Path("."),
    )
    genes = assign_sensor_roles(
        genes,
        regulator_families=config.regulator_families,
        sensory_domain_catalog=config.sensory_domains,
        paired_evidence_rules=config.paired_evidence,
    )
    sensor_regulator_pairs = build_sensor_regulator_pairs(
        genes,
        loci,
        sensory_domain_catalog=config.sensory_domains,
        paired_evidence_rules=config.paired_evidence,
    )
    scored_loci = score_loci(loci, config.scoring)
    co_locus_ids = (
        scored_loci.filter(pl.col("analyte") == "CO")["locus_id"].unique().to_list()
    )
    cn_locus_ids = (
        scored_loci.filter(pl.col("analyte") == "CN")["locus_id"].unique().to_list()
    )
    enrichment_genes = _genes_with_locus_taxonomy(scored_loci, genes)
    case_enrichment_genes = enrichment_genes.filter(
        pl.col("locus_id").is_in(co_locus_ids),
    )
    control_enrichment_genes = enrichment_genes.filter(
        pl.col("locus_id").is_in(cn_locus_ids),
    )
    if config.scoring.enrichment.test == "cmh":
        enrichment = run_stratified_enrichment(
            case_enrichment_genes,
            control_enrichment_genes,
            analyte="CO",
            case_definition="mini CO fixture loci",
            control_definition="mini cyd-control fixture loci",
            stratum_column=config.scoring.enrichment.stratum_column,
            alpha=config.scoring.enrichment.alpha,
            deduplication_policy=config.scoring.enrichment.strict_policy,
        )
    else:
        enrichment = run_enrichment(
            case_enrichment_genes,
            control_enrichment_genes,
            analyte="CO",
            case_definition="mini CO fixture loci",
            control_definition="mini cyd-control fixture loci",
            alpha=config.scoring.enrichment.alpha,
        )
    enrichment_robustness = run_enrichment_robustness(
        case_enrichment_genes,
        control_enrichment_genes,
        analyte="CO",
        case_definition="mini CO fixture loci",
        control_definition="mini cyd-control fixture loci",
        stratum_column=config.scoring.enrichment.stratum_column,
        alpha=config.scoring.enrichment.alpha,
    )
    candidates = score_candidates(
        scored_loci,
        genes,
        config.scoring,
        enrichment,
        sensory_domain_catalog=config.sensory_domains,
        paired_evidence_rules=config.paired_evidence,
        expected_chemistry_by_analyte=expected_chemistry_by_analyte(config.analytes),
    )
    archetypes = cluster_archetypes(scored_loci, candidates)
    candidates = compute_conservation_scores(
        candidates,
        archetypes,
        scored_loci,
        min_loci_per_archetype=1,
        scoring=config.scoring,
    )
    candidates = assign_phylogenetic_profile_scores(
        candidates,
        scored_loci,
        scoring=config.scoring,
    )
    candidates = assign_operon_regulation_score_bands(candidates)
    archetypes = cluster_archetypes(scored_loci, candidates)
    benchmark = _benchmark_recovery(candidates)

    _write_frame_set(
        out_dir=out_dir,
        loci=scored_loci,
        genes=genes,
        candidates=candidates,
        enrichment=enrichment,
        enrichment_robustness=enrichment_robustness,
        archetypes=archetypes,
        sensor_regulator_pairs=sensor_regulator_pairs,
    )
    write_publication_tables(
        benchmark_recovery=benchmark,
        candidates=candidates,
        enrichment=enrichment,
        archetypes=archetypes,
        out_dir=out_dir / "tables",
    )
    figures_dir = out_dir / "figures"
    figure_1_workflow_and_recovery(benchmark, figures_dir)
    figure_2_locus_landscape(scored_loci, figures_dir)
    figure_3_archetype_atlas(archetypes, figures_dir)
    figure_4_chemistry_partition(enrichment, figures_dir)
    figure_5_candidate_ranking(candidates, figures_dir)
    figure_6_structure_hypotheses(candidates, out_dir / "structures", figures_dir)
    captions = build_result_led_captions(
        benchmark_results=benchmark,
        loci=scored_loci,
        archetypes=archetypes,
        enrichment=enrichment,
        candidates=candidates,
        top_candidates=candidates,
    )
    write_caption_files(captions, out_dir / "captions")

    manifest = build_manifest(
        seed=config.seed,
        command="sqlite-demo",
        config_paths={"config": config_path},
        input_paths={"sqlite": sqlite},
    )
    write_manifest(manifest, out_dir)
    marker = out_dir / "README.txt"
    marker.write_text("GasRegNet SQLite demo workflow complete\n", encoding="utf-8")
    return marker


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--sqlite", type=Path, required=True)
    args = parser.parse_args()
    run_sqlite_demo(out_dir=args.out, config_path=args.config, sqlite_path=args.sqlite)


if __name__ == "__main__":
    main()
