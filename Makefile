.PHONY: sync lint test repro repro-real clean

sync:
	uv sync --extra dev

lint:
	uv run ruff check gasregnet
	uv run mypy --strict gasregnet

test:
	uv run pytest -q --cov=gasregnet

repro:
	uv run snakemake -s workflows/sqlite_mode.smk --cores 1

repro-real:
	uv run snakemake -s workflows/full_discovery.smk --cores 1

clean:
	rm -rf results/* cache .snakemake .pytest_cache .ruff_cache .mypy_cache
	touch results/.gitkeep
