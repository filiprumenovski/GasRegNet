#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORK_DIR="${GASREGNET_ACCEPTANCE_WORKDIR:-${ROOT_DIR}/results/scaling_probe_dryrun}"
CORES="${CORES:-1}"

cd "$ROOT_DIR"
rm -rf "$WORK_DIR"
uv run python scripts/generate_scaling_corpus.py \
  --out-dir "$WORK_DIR" \
  --datasets "${GASREGNET_ACCEPTANCE_DATASETS:-2}" \
  --genes-per-dataset "${GASREGNET_ACCEPTANCE_GENES_PER_DATASET:-24}" \
  --window-genes "${GASREGNET_ACCEPTANCE_WINDOW_GENES:-4}"

uv run snakemake \
  --snakefile workflows/corpus_discovery.smk \
  --dry-run \
  --cores "$CORES" \
  --config \
  out_dir="$WORK_DIR/results/corpus" \
  corpus_config="$WORK_DIR/configs/corpus_discovery.yaml" \
  config_path="configs/headline.yaml" \
  root="$WORK_DIR" \
  catalogs="$WORK_DIR/configs/refseq_catalogs.yaml" \
  scan_config="$WORK_DIR/configs/refseq_scan.yaml" \
  mode="smoke"
