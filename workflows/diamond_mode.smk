OUT_DIR = config.get("out_dir", "results/diamond")
QUERY = config.get("query", "data/seeds/co_anchor_seeds.faa")
DB = config.get("db", "databases/bacteria.dmnd")


rule all:
    input:
        f"{OUT_DIR}/hits.parquet",


rule diamond_search:
    input:
        query=QUERY,
        db=DB,
    output:
        hits=f"{OUT_DIR}/hits.parquet",
    shell:
        """
        uv run gasregnet diamond-search \
          --query {input.query} \
          --db {input.db} \
          --out {output.hits}
        """
