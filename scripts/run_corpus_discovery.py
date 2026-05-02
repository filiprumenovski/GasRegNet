"""Run the RefSeq corpus discovery pipeline."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import polars as pl
import yaml  # type: ignore[import-untyped]

from gasregnet.annotation.ecology import score_taxonomic_context_by_analyte
from gasregnet.annotation.regulators import classify_regulators
from gasregnet.annotation.roles import (
    assign_sensor_roles,
    build_sensor_regulator_pairs,
)
from gasregnet.archetypes.cluster import cluster_archetypes
from gasregnet.benchmark import evaluate_benchmark
from gasregnet.config import load_config, resolve_and_dump
from gasregnet.datasets.refseq import (
    detect_refseq_anchor_hits,
    extract_refseq_neighborhoods,
)
from gasregnet.io.parquet import write_parquet
from gasregnet.manifest import build_manifest, write_manifest
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
    AnchorHitsSchema,
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


def _read_corpus_config(path: Path) -> dict[str, Any]:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"corpus config must be a mapping: {path}")
    return payload


def _write_frame_set(
    *,
    out_dir: Path,
    anchor_hits: pl.DataFrame,
    loci: pl.DataFrame,
    genes: pl.DataFrame,
    candidates: pl.DataFrame,
    enrichment: pl.DataFrame,
    enrichment_robustness: pl.DataFrame,
    archetypes: pl.DataFrame,
    sensor_regulator_pairs: pl.DataFrame,
) -> None:
    intermediate = out_dir / "intermediate"
    write_parquet(anchor_hits, intermediate / "anchor_hits.parquet", AnchorHitsSchema)
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


def _enrichment(
    loci: pl.DataFrame,
    genes: pl.DataFrame,
    *,
    test: str,
    alpha: float,
    stratum_column: str,
    strict_policy: str,
) -> pl.DataFrame:
    genes = _genes_with_locus_taxonomy(loci, genes)
    case_loci = loci.filter(pl.col("analyte") == "CO")["locus_id"].to_list()
    control_loci = loci.filter(pl.col("analyte") == "CN")["locus_id"].to_list()
    case_genes = genes.filter(pl.col("locus_id").is_in(case_loci))
    control_genes = genes.filter(pl.col("locus_id").is_in(control_loci))
    if test == "cmh":
        return run_stratified_enrichment(
            case_genes,
            control_genes,
            analyte="CO",
            case_definition="RefSeq corpus CO anchor neighborhoods",
            control_definition="RefSeq corpus cyd-control neighborhoods",
            stratum_column=stratum_column,
            alpha=alpha,
            deduplication_policy=strict_policy,
        )
    return run_enrichment(
        case_genes,
        control_genes,
        analyte="CO",
        case_definition="RefSeq corpus CO anchor neighborhoods",
        control_definition="RefSeq corpus cyd-control neighborhoods",
        alpha=alpha,
    )


def _enrichment_robustness(
    loci: pl.DataFrame,
    genes: pl.DataFrame,
    *,
    alpha: float,
    stratum_column: str,
) -> pl.DataFrame:
    genes = _genes_with_locus_taxonomy(loci, genes)
    case_loci = loci.filter(pl.col("analyte") == "CO")["locus_id"].to_list()
    control_loci = loci.filter(pl.col("analyte") == "CN")["locus_id"].to_list()
    return run_enrichment_robustness(
        genes.filter(pl.col("locus_id").is_in(case_loci)),
        genes.filter(pl.col("locus_id").is_in(control_loci)),
        analyte="CO",
        case_definition="RefSeq corpus CO anchor neighborhoods",
        control_definition="RefSeq corpus cyd-control neighborhoods",
        stratum_column=stratum_column,
        alpha=alpha,
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


def run_corpus_discovery(
    *,
    out_dir: Path,
    corpus_config_path: Path,
    config_path: Path,
    root: Path,
) -> Path:
    """Run corpus-scale discovery from RefSeq catalogs to report artifacts."""

    ensure_out_dir(out_dir)
    corpus_config = _read_corpus_config(corpus_config_path)
    config = load_config(config_path)
    resolve_and_dump(config, out_dir)

    catalogs = root / str(corpus_config["catalogs"])
    scan_config = root / str(corpus_config["scan_config"])
    benchmark_csv = root / str(corpus_config["benchmark"])
    window_genes = int(corpus_config.get("window_genes", 10))
    mode = str(corpus_config.get("mode", "smoke"))
    profile_dir = root / str(corpus_config.get("profile_dir", "data/profiles"))
    bitscore_threshold = corpus_config.get("bitscore_threshold")
    e_value_threshold = float(corpus_config.get("e_value_threshold", 1e-20))

    anchor_hits = detect_refseq_anchor_hits(
        catalogs,
        scan_config,
        root=root,
        mode=mode,
        config=config_path,
        profile_dir=profile_dir,
        bitscore_threshold=(
            None if bitscore_threshold is None else float(bitscore_threshold)
        ),
        e_value_threshold=e_value_threshold,
    )
    loci, genes = extract_refseq_neighborhoods(
        anchor_hits,
        catalogs,
        root=root,
        window_genes=window_genes,
    )
    loci = score_taxonomic_context_by_analyte(loci, config.analytes, root=root)
    genes = classify_regulators(genes, config.regulator_families)
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
    enrichment = _enrichment(
        scored_loci,
        genes,
        test=config.scoring.enrichment.test,
        alpha=config.scoring.enrichment.alpha,
        stratum_column=config.scoring.enrichment.stratum_column,
        strict_policy=config.scoring.enrichment.strict_policy,
    )
    enrichment_robustness = _enrichment_robustness(
        scored_loci,
        genes,
        alpha=config.scoring.enrichment.alpha,
        stratum_column=config.scoring.enrichment.stratum_column,
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
        scoring=config.scoring,
    )
    candidates = assign_phylogenetic_profile_scores(
        candidates,
        scored_loci,
        scoring=config.scoring,
    )
    candidates = assign_operon_regulation_score_bands(candidates)
    archetypes = cluster_archetypes(scored_loci, candidates)
    benchmark = evaluate_benchmark(benchmark_csv, anchor_hits, candidates)

    _write_frame_set(
        out_dir=out_dir,
        anchor_hits=anchor_hits,
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
    write_manifest(
        build_manifest(
            seed=config.seed,
            command="corpus-discovery",
            config_paths={
                "corpus_config": corpus_config_path,
                "config": config_path,
                "catalogs": catalogs,
                "scan_config": scan_config,
            },
            input_paths={"benchmark": benchmark_csv},
        ),
        out_dir,
    )
    marker = out_dir / "README.txt"
    marker.write_text(
        "GasRegNet RefSeq corpus discovery workflow complete\n",
        encoding="utf-8",
    )
    return marker


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--corpus-config", type=Path, required=True)
    parser.add_argument("--config", type=Path, default=Path("configs/headline.yaml"))
    parser.add_argument("--root", type=Path, default=Path("."))
    args = parser.parse_args()
    run_corpus_discovery(
        out_dir=args.out,
        corpus_config_path=args.corpus_config,
        config_path=args.config,
        root=args.root,
    )


if __name__ == "__main__":
    main()
