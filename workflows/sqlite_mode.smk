configfile: "configs/headline.yaml"

OUT_DIR = config.get("out_dir", "results/repro")
CONFIG_PATH = config.get("config_path", "configs/headline.yaml")
SQLITE = config.get("sqlite", "tests/fixtures/mini_efi.sqlite")


rule all:
    input:
        f"{OUT_DIR}/README.txt",
        f"{OUT_DIR}/manifest.json",
        f"{OUT_DIR}/config.resolved.yaml",
        f"{OUT_DIR}/intermediate/loci.parquet",
        f"{OUT_DIR}/intermediate/genes.parquet",
        f"{OUT_DIR}/intermediate/sensor_regulator_pairs.parquet",
        f"{OUT_DIR}/intermediate/candidates.parquet",
        f"{OUT_DIR}/intermediate/enrichment.parquet",
        f"{OUT_DIR}/intermediate/archetypes.parquet",
        f"{OUT_DIR}/tables/T1_benchmark_recovery.csv",
        f"{OUT_DIR}/tables/T6_tool_feature_comparison.md",
        f"{OUT_DIR}/figures/figure_1_workflow_and_recovery.png",
        f"{OUT_DIR}/figures/figure_6_structure_hypotheses.svg",
        f"{OUT_DIR}/captions/figure_4_chemistry_partition.md",


rule sqlite_demo:
    output:
        marker=f"{OUT_DIR}/README.txt",
        manifest=f"{OUT_DIR}/manifest.json",
        resolved_config=f"{OUT_DIR}/config.resolved.yaml",
        loci=f"{OUT_DIR}/intermediate/loci.parquet",
        genes=f"{OUT_DIR}/intermediate/genes.parquet",
        pairs=f"{OUT_DIR}/intermediate/sensor_regulator_pairs.parquet",
        candidates=f"{OUT_DIR}/intermediate/candidates.parquet",
        enrichment=f"{OUT_DIR}/intermediate/enrichment.parquet",
        archetypes=f"{OUT_DIR}/intermediate/archetypes.parquet",
        benchmark_table=f"{OUT_DIR}/tables/T1_benchmark_recovery.csv",
        comparison_table=f"{OUT_DIR}/tables/T6_tool_feature_comparison.md",
        recovery_figure=f"{OUT_DIR}/figures/figure_1_workflow_and_recovery.png",
        structure_figure=f"{OUT_DIR}/figures/figure_6_structure_hypotheses.svg",
        partition_caption=f"{OUT_DIR}/captions/figure_4_chemistry_partition.md",
    params:
        out_dir=OUT_DIR,
        config=CONFIG_PATH,
        sqlite=SQLITE,
    shell:
        """
        uv run python scripts/run_sqlite_demo.py \
          --out {params.out_dir} \
          --config {params.config} \
          --sqlite {params.sqlite}
        """
