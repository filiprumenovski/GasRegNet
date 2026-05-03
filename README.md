# GasRegNet

GasRegNet is a comparative genomics workflow for finding candidate bacterial
gas-response regulators from genome neighborhoods. It starts from functional
anchor genes, retrieves nearby genes, annotates regulator and sensor evidence,
scores candidate regulators, and writes reproducible Parquet, table, figure, and
manifest outputs.

The current analyte set covers carbon monoxide (`CO`), nitric oxide (`NO`),
oxygen (`O2`), cyanide handling (`CN`), and a separate cytochrome bd respiration
control (`cyd_control`). The control split is intentional: `cydA/cydB/cydX`
neighborhoods are not reported as cyanide-sensing headline biology.

## What It Does

- Imports EFI-GNT SQLite fixtures and RefSeq FASTA/GFF catalogs.
- Builds and validates typed Polars/Pandera data products.
- Detects anchor proteins with profile HMMs, DIAMOND seed rescue, and explicit
  smoke-mode term scans.
- Extracts genome-resolved neighborhoods around anchors.
- Classifies nearby regulator and sensor candidates.
- Scores candidates with decomposable evidence components.
- Runs matched-control enrichment tests.
- Clusters recurrent locus-neighborhood archetypes.
- Writes reproducible reports, manifests, and readiness checks.

## What It Does Not Claim

- `candidate_score` is a deterministic prioritization score, not a calibrated
  probability that a protein regulates gas response.
- The score-band output is explicitly uncalibrated.
- Archetyping is locus-neighborhood archetyping, not protein/fold archetyping.
- UniRef90 Parquets are useful for sequence-scale anchor discovery, but they do
  not contain genome coordinates and cannot replace RefSeq/GTDB-style genome
  data for neighborhood/regulator extraction.

## Install

GasRegNet targets Python `3.11.9` and uses `uv`.

```bash
uv sync --extra dev
```

External binaries are only needed for specific paths:

- HMMER: `hmmsearch`
- DIAMOND: `diamond`
- Optional: MAFFT, Foldseek, MEME-suite tools

Check local tool resolution:

```bash
make check-tools
```

## Quickstart

Run the deterministic SQLite fixture workflow:

```bash
make repro
```

Run tests and type checks:

```bash
make lint
uv run pytest -q
```

Run the UniRef90 readiness gate:

```bash
make uniref90-readiness
```

This verifies tool resolution, builds per-family DIAMOND seed databases, runs
targeted schema/profile/store tests, and dry-runs the sharded corpus workflow.

## RefSeq Corpus Workflow

Fetch configured external assets:

```bash
make assets
make datasets
```

Index configured RefSeq catalogs:

```bash
make index-datasets
make summarize-datasets
```

Run the small corpus workflow:

```bash
make corpus-repro
```

Corpus outputs are written under `results/corpus/`. Runtime data, external
downloads, generated DIAMOND databases, corpus stores, and DuckDB catalogs are
ignored by Git.

## Main CLI Commands

```bash
uv run gasregnet validate-config --config configs --out results/config_check
uv run gasregnet build-profiles --config configs --out-dir data/profiles
uv run gasregnet build-seed-databases --config configs/headline.yaml --out data/profiles/diamond_seeds
uv run gasregnet index-refseq-corpus-store --manifest configs/refseq_catalogs.yaml --store data/corpus_store
uv run gasregnet enumerate-shards --store data/corpus_store --out data/corpus_store/shards.parquet
uv run gasregnet detect-anchors-shard --store data/corpus_store --shards data/corpus_store/shards.parquet --shard-id <shard> --seed-diamond-dir data/profiles/diamond_seeds --out results/corpus/intermediate/anchor_hits.<shard>.parquet
uv run gasregnet extract-neighborhoods-shard --store data/corpus_store --anchor-hits <anchor_hits.parquet> --out results/corpus
uv run gasregnet report --results results/corpus --out results/corpus
```

Use `uv run gasregnet --help` for the full command list.

## Data Products

Canonical intermediate tables are Parquet:

- `anchor_hits.parquet`
- `loci.parquet`
- `genes.parquet`
- `sensor_regulator_pairs.parquet`
- `candidates.parquet`
- `enrichment.parquet`
- `archetypes.parquet`

Report outputs include CSV/Markdown tables, PNG/SVG figures, captions, resolved
config, and a run manifest.

## Development Checks

```bash
uv run ruff check gasregnet scripts tests
uv run mypy --strict gasregnet scripts/fetch_assets.py scripts/build_test_fixtures.py scripts/build_profiles.py scripts/run_corpus_discovery.py
uv run pytest -q
./scripts/run_acceptance.sh
```

The repository is intended to stay clean after these commands: generated
outputs belong under ignored runtime directories, not in source control.
