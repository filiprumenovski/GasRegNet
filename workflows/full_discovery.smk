configfile: "configs/headline.yaml"

OUT_DIR = config.get("out_dir", "results/full_discovery")
CONFIG_PATH = config.get("config_path", "configs/headline.yaml")
SQLITE = config.get("sqlite", "tests/fixtures/mini_efi.sqlite")


rule all:
    input:
        f"{OUT_DIR}/README.txt",
        f"{OUT_DIR}/manifest.json",
        f"{OUT_DIR}/intermediate/candidates.parquet",
        f"{OUT_DIR}/tables/T6_tool_feature_comparison.csv",
        f"{OUT_DIR}/figures/figure_4_chemistry_partition.png",
        f"{OUT_DIR}/captions/figure_5_candidate_ranking.md",


rule fixture_full_discovery:
    output:
        marker=f"{OUT_DIR}/README.txt",
        manifest=f"{OUT_DIR}/manifest.json",
        candidates=f"{OUT_DIR}/intermediate/candidates.parquet",
        comparison_table=f"{OUT_DIR}/tables/T6_tool_feature_comparison.csv",
        partition_figure=f"{OUT_DIR}/figures/figure_4_chemistry_partition.png",
        candidate_caption=f"{OUT_DIR}/captions/figure_5_candidate_ranking.md",
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
