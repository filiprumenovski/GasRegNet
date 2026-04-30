configfile: "configs/corpus_discovery.yaml"

OUT_DIR = config.get("out_dir", "results/corpus")
CORPUS_CONFIG = config.get("corpus_config", "configs/corpus_discovery.yaml")
HEADLINE_CONFIG = config.get("config_path", "configs/headline.yaml")
ROOT = config.get("root", ".")


rule all:
    input:
        f"{OUT_DIR}/README.txt",
        f"{OUT_DIR}/manifest.json",
        f"{OUT_DIR}/intermediate/anchor_hits.parquet",
        f"{OUT_DIR}/intermediate/loci.parquet",
        f"{OUT_DIR}/intermediate/genes.parquet",
        f"{OUT_DIR}/intermediate/sensor_regulator_pairs.parquet",
        f"{OUT_DIR}/intermediate/candidates.parquet",
        f"{OUT_DIR}/intermediate/enrichment.parquet",
        f"{OUT_DIR}/intermediate/enrichment_robustness.parquet",
        f"{OUT_DIR}/tables/T1_benchmark_recovery.csv",
        f"{OUT_DIR}/tables/T6_tool_feature_comparison.csv",
        f"{OUT_DIR}/figures/figure_1_workflow_and_recovery.png",
        f"{OUT_DIR}/captions/figure_5_candidate_ranking.md",


rule corpus_discovery:
    output:
        marker=f"{OUT_DIR}/README.txt",
        manifest=f"{OUT_DIR}/manifest.json",
        anchors=f"{OUT_DIR}/intermediate/anchor_hits.parquet",
        loci=f"{OUT_DIR}/intermediate/loci.parquet",
        genes=f"{OUT_DIR}/intermediate/genes.parquet",
        pairs=f"{OUT_DIR}/intermediate/sensor_regulator_pairs.parquet",
        candidates=f"{OUT_DIR}/intermediate/candidates.parquet",
        enrichment=f"{OUT_DIR}/intermediate/enrichment.parquet",
        enrichment_robustness=f"{OUT_DIR}/intermediate/enrichment_robustness.parquet",
        benchmark=f"{OUT_DIR}/tables/T1_benchmark_recovery.csv",
        comparison_table=f"{OUT_DIR}/tables/T6_tool_feature_comparison.csv",
        recovery_figure=f"{OUT_DIR}/figures/figure_1_workflow_and_recovery.png",
        candidate_caption=f"{OUT_DIR}/captions/figure_5_candidate_ranking.md",
    params:
        out_dir=OUT_DIR,
        corpus_config=CORPUS_CONFIG,
        headline_config=HEADLINE_CONFIG,
        root=ROOT,
    shell:
        """
        uv run python scripts/run_corpus_discovery.py \
          --out {params.out_dir} \
          --corpus-config {params.corpus_config} \
          --config {params.headline_config} \
          --root {params.root}
        """
