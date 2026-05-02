configfile: "configs/corpus_discovery.yaml"

OUT_DIR = config.get("out_dir", "results/corpus")
CORPUS_CONFIG = config.get("corpus_config", "configs/corpus_discovery.yaml")
HEADLINE_CONFIG = config.get("config_path", "configs/headline.yaml")
ROOT = config.get("root", ".")
STORE = config.get("store", "data/corpus_store")
SEED_DIAMOND_DIR = config.get("seed_diamond_dir", "")
SEED_DIAMOND_ARG = "" if not SEED_DIAMOND_DIR else f"--seed-diamond-dir {SEED_DIAMOND_DIR}"
SHARD_CONFIG = config.get("sharding", {})
SHARD_CONFIG = SHARD_CONFIG if isinstance(SHARD_CONFIG, dict) else {}
SHARD_STRATEGY = config.get("sharding_strategy", SHARD_CONFIG.get("strategy", "by_phylum"))
N_GENOMES_PER_SHARD = int(
    config.get(
        "n_genomes_per_shard",
        SHARD_CONFIG.get("n_genomes_per_shard", 1000),
    ),
)
MAX_SHARDS = config.get("max_shards", SHARD_CONFIG.get("max_shards"))
MAX_SHARDS = None if MAX_SHARDS in (None, "null", "") else int(MAX_SHARDS)
SHARDS_PATH = config.get("shards", f"{STORE}/shards.parquet")


def resolve_shards():
    from pathlib import Path

    from gasregnet.datasets.sharding import read_shards, write_shards

    shards_path = Path(SHARDS_PATH)
    store = Path(STORE)
    if (store / "datasets.parquet").exists():
        write_shards(
            store,
            shards_path,
            SHARD_STRATEGY,
            n_genomes_per_shard=N_GENOMES_PER_SHARD,
            max_shards=MAX_SHARDS,
        )
    if shards_path.exists():
        return [shard.shard_id for shard in read_shards(shards_path)]
    return SHARD_CONFIG.get("shards", ["all"])


SHARDS = resolve_shards()


rule all:
    input:
        SHARDS_PATH,
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
        f"{OUT_DIR}/intermediate/sharded/anchor_hits.parquet",
        f"{OUT_DIR}/intermediate/sharded/loci.parquet",
        f"{OUT_DIR}/intermediate/sharded/genes.parquet",


rule enumerate_shards:
    output:
        SHARDS_PATH,
    params:
        store=STORE,
        strategy=SHARD_STRATEGY,
        n_genomes_per_shard=N_GENOMES_PER_SHARD,
        max_shards_arg="" if MAX_SHARDS is None else f"--max-shards {MAX_SHARDS}",
    shell:
        """
        uv run gasregnet enumerate-shards \
          --store {params.store} \
          --strategy {params.strategy} \
          --n-genomes-per-shard {params.n_genomes_per_shard} \
          {params.max_shards_arg} \
          --out {output}
        """


rule detect_anchors_shard:
    input:
        shards=SHARDS_PATH,
    output:
        f"{OUT_DIR}/intermediate/sharded/anchor_hits.{{shard}}.parquet",
    params:
        config_path=HEADLINE_CONFIG,
        profile_dir=config.get("profile_dir", "data/profiles"),
        e_value_threshold=config.get("e_value_threshold", 1e-20),
        store=STORE,
        shards=SHARDS_PATH,
        seed_diamond_arg=SEED_DIAMOND_ARG,
    shell:
        """
        uv run gasregnet detect-anchors-shard \
          --store {params.store} \
          --shards {params.shards} \
          --shard-id {wildcards.shard} \
          --config {params.config_path} \
          --profile-dir {params.profile_dir} \
          --e-value-threshold {params.e_value_threshold} \
          {params.seed_diamond_arg} \
          --out {output}
        """


rule gather_anchor_shards:
    input:
        expand(f"{OUT_DIR}/intermediate/sharded/anchor_hits.{{shard}}.parquet", shard=SHARDS),
    output:
        f"{OUT_DIR}/intermediate/sharded/anchor_hits.parquet",
    run:
        import polars as pl

        frames = [pl.read_parquet(path) for path in input]
        pl.concat(frames, how="vertical").unique().write_parquet(output[0])


rule extract_neighborhoods_shard:
    input:
        f"{OUT_DIR}/intermediate/sharded/anchor_hits.{{shard}}.parquet",
    output:
        loci=f"{OUT_DIR}/intermediate/sharded/loci.{{shard}}.parquet",
        genes=f"{OUT_DIR}/intermediate/sharded/genes.{{shard}}.parquet",
    params:
        store=STORE,
        window_genes=config.get("window_genes", 10),
        out_dir=f"{OUT_DIR}/intermediate/sharded/extract.{{shard}}",
    shell:
        """
        uv run gasregnet extract-neighborhoods-shard \
          --anchor-hits {input} \
          --store {params.store} \
          --window-genes {params.window_genes} \
          --out {params.out_dir}
        cp {params.out_dir}/intermediate/loci.parquet {output.loci}
        cp {params.out_dir}/intermediate/genes.parquet {output.genes}
        """


rule gather_neighborhood_shards:
    input:
        loci=expand(f"{OUT_DIR}/intermediate/sharded/loci.{{shard}}.parquet", shard=SHARDS),
        genes=expand(f"{OUT_DIR}/intermediate/sharded/genes.{{shard}}.parquet", shard=SHARDS),
    output:
        loci=f"{OUT_DIR}/intermediate/sharded/loci.parquet",
        genes=f"{OUT_DIR}/intermediate/sharded/genes.parquet",
    run:
        import polars as pl

        pl.concat([pl.read_parquet(path) for path in input.loci], how="vertical").write_parquet(output.loci)
        pl.concat([pl.read_parquet(path) for path in input.genes], how="vertical").write_parquet(output.genes)


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
