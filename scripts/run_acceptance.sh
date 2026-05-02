#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

export GASREGNET_ACCEPTANCE="${GASREGNET_ACCEPTANCE:-1}"
uv run pytest tests/integration/scaling_probe "$@"
