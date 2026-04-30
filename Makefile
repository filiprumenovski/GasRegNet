.PHONY: sync lint test assets repro repro-real clean

sync:
	uv sync --extra dev

lint:
	uv run ruff check gasregnet scripts tests
	uv run mypy --strict gasregnet scripts/fetch_assets.py scripts/build_test_fixtures.py

test:
	uv run pytest -q --cov=gasregnet

assets:
	uv run python scripts/fetch_assets.py --manifest configs/assets.yaml --force

repro:
	uv run snakemake -s workflows/sqlite_mode.smk --cores 1

repro-real:
	uv run snakemake -s workflows/full_discovery.smk --cores 1

clean:
	rm -rf results/* cache .snakemake .pytest_cache .ruff_cache .mypy_cache
	touch results/.gitkeep
