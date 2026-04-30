# GasRegNet

GasRegNet is a comparative genomics package for discovering candidate bacterial
gas-sensing transcriptional regulators. It implements the project contracts in
`ARCHITECTURE.md`, `PREPRINT_PLAN.md`, and `ENGINEERING_PROMPT.md`: validated
configuration, EFI-GNT SQLite ingestion, annotation, scoring, enrichment,
archetype clustering, report tables, figures, captions, manifests, and a local
Snakemake reproducibility path.

## Quickstart

```bash
make sync
make assets
make datasets
make index-datasets
make lint
make test
make repro
```

`make repro` runs the deterministic mini EFI-GNT SQLite fixture through the
implemented pipeline and writes a complete local report bundle under
`results/repro/`:

- `intermediate/*.parquet`
- `tables/T1_*.csv` through `tables/T6_*.csv` plus Markdown copies
- six `figures/*.png` and `figures/*.svg` files
- six result-led `captions/*.md` files
- `manifest.json` and `config.resolved.yaml`

The fixture demonstrates the engineering path. Publication-scale biological
claims still require the curated benchmark and real EFI-GNT input assets.

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

## CLI

```bash
uv run gasregnet validate-config --config configs --out results/config_check
uv run gasregnet build-benchmark --out data/benchmarks/benchmark_v1.csv
uv run gasregnet run-sqlite --sqlite tests/fixtures/mini_efi.sqlite --analytes CO --analytes CN --config configs --out results/sqlite_demo
uv run gasregnet diamond-search --query data/seeds/co_anchor_seeds.faa --db databases/bacteria.dmnd --out cache/co_diamond_hits.parquet
```

External search commands require their corresponding local databases and
binaries. The checked-in reproducibility workflow does not require network
access or external biological databases.
