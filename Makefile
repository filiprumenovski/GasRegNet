.PHONY: sync lint test assets datasets index-datasets summarize-datasets scan-datasets repro repro-real clean

sync:
	uv sync --extra dev

lint:
	uv run ruff check gasregnet scripts tests
	uv run mypy --strict gasregnet scripts/fetch_assets.py scripts/build_test_fixtures.py

test:
	uv run pytest -q --cov=gasregnet

assets:
	uv run gasregnet fetch-assets --manifest configs/assets.yaml --force

datasets:
	uv run gasregnet fetch-assets --manifest configs/datasets.yaml --downloader aria2 --force

index-datasets:
	uv run gasregnet index-refseq-corpus --manifest configs/refseq_catalogs.yaml

summarize-datasets:
	uv run gasregnet summarize-refseq-corpus --manifest configs/refseq_catalogs.yaml

scan-datasets:
	uv run gasregnet scan-refseq-corpus \
		--manifest configs/refseq_catalogs.yaml \
		--scan-config configs/refseq_scan.yaml \
		--out results/refseq_anchor_scan.csv

repro:
	uv run snakemake -s workflows/sqlite_mode.smk --cores 1

repro-real:
	uv run snakemake -s workflows/full_discovery.smk --cores 1

clean:
	rm -rf results/* cache .snakemake .pytest_cache .ruff_cache .mypy_cache
	touch results/.gitkeep
