# GasRegNet

GasRegNet is a comparative genomics package for discovering candidate bacterial
gas-sensing transcriptional regulators across CO, NO, and O2, with cyd loci
kept as a cyanide-resistant respiration control rather than a cyanide-sensing
headline. It implements the project contracts in
`ARCHITECTURE.md`, `PREPRINT_PLAN.md`, and `ENGINEERING_PROMPT.md`: validated
configuration, EFI-GNT SQLite ingestion, DuckDB RefSeq catalogs, corpus anchor
detection, RefSeq neighborhood extraction, annotation, scoring, deterministic
operon-level regulation score bands, phylogenetic profile
co-occurrence evidence, synthetic-truth calibration fixtures, enrichment,
archetype clustering, report tables, figures, captions, manifests, and local
Snakemake reproducibility paths.

## Quickstart

```bash
make sync
make assets
make datasets
make index-datasets
make lint
make test
make repro
make corpus-repro
```

`make repro` runs the deterministic mini EFI-GNT SQLite fixture through the
implemented pipeline and writes a complete local report bundle under
`results/repro/`:

- `intermediate/*.parquet`
- `tables/T1_*.csv` through `tables/T6_*.csv` plus Markdown copies
- six `figures/*.png` and `figures/*.svg` files
- six result-led `captions/*.md` files
- `manifest.json` and `config.resolved.yaml`

The fixture demonstrates the engineering path. `make corpus-repro` downloads the
configured RefSeq assets, indexes them into DuckDB, detects smoke-mode anchors,
extracts RefSeq neighborhoods, and writes a report bundle under
`results/corpus/`. Profile-mode anchor detection uses reproducible HMM profiles
with seed back-confirmation where local seed evidence supports the hit.

## Development

The package targets Python 3.11.9 and uses `uv` for environment management.

```bash
uv sync
uv run ruff check gasregnet
uv run mypy --strict gasregnet
uv run pytest -q --cov=gasregnet
```

Runtime outputs are written under `results/` and are ignored by version control.
External seed FASTAs are declared in `configs/assets.yaml` and can be refreshed
with `make assets`.
Larger local datasets are declared in `configs/datasets.yaml`, downloaded with
`aria2c` by `make datasets`, and written under ignored `data/external/`.
`make index-datasets` converts those assets into local DuckDB reference catalogs
under ignored `databases/`.
`make summarize-datasets` reports catalog sizes and feature-linkage counts.
`make scan-datasets` scans indexed catalogs for configured gas-anchor terms and
writes `results/refseq_anchor_scan.csv`.
`make corpus-repro` runs the small RefSeq corpus workflow end to end and writes
`anchor_hits.parquet`, canonical `loci.parquet` and `genes.parquet`, scored
tables, candidate score bands, figures, captions, and a run
manifest under `results/corpus/`.

## CLI

```bash
uv run gasregnet validate-config --config configs --out results/config_check
uv run gasregnet build-benchmark --version v2 --out data/benchmarks/regulators_v2.csv
uv run gasregnet detect-anchors --out results/corpus/intermediate/anchor_hits.parquet
uv run gasregnet extract-neighborhoods --anchor-hits results/corpus/intermediate/anchor_hits.parquet --out results/corpus
uv run gasregnet evaluate-benchmark --anchor-hits results/corpus/intermediate/anchor_hits.parquet --out results/corpus/tables/T1_benchmark_recovery.csv
uv run gasregnet simulate-synthetic-truth --out results/synthetic_truth --n-genomes 24 --annotation-noise 0.1
uv run gasregnet run-sqlite --sqlite tests/fixtures/mini_efi.sqlite --analytes CO --analytes CN --config configs --out results/sqlite_demo
uv run gasregnet diamond-search --query data/seeds/co_anchor_seeds.faa --db databases/bacteria.dmnd --out cache/co_diamond_hits.parquet
uv run gasregnet corpus-discovery --out results/corpus
```

External search commands require their corresponding local databases and
binaries. The checked-in reproducibility workflow does not require network
access or external biological databases.
