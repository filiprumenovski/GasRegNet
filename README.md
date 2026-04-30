# GasRegNet

GasRegNet is a comparative genomics package for discovering candidate bacterial gas-sensing transcriptional regulators. It is being scaffolded from the project contracts in `ARCHITECTURE.md`, `PREPRINT_PLAN.md`, and `ENGINEERING_PROMPT.md`.

## Quickstart

```bash
make sync
make lint
make test
make repro
```

The initial scaffold contains the repository layout, pinned Python project metadata, workflow entry points, configuration placeholders, and empty module boundaries. Implementation follows the task graph in `ENGINEERING_PROMPT.md`.

## Development

The package targets Python 3.11.9 and uses `uv` for environment management.

```bash
uv sync
uv run ruff check gasregnet
uv run mypy --strict gasregnet
uv run pytest -q --cov=gasregnet
```

Runtime outputs are written under `results/` and are ignored by version control.
