.PHONY: sync lint test assets profiles seed-databases datasets index-datasets summarize-datasets scan-datasets discover-motifs embed-proteins foldseek-search fit-bayesian build-views emit-provenance check-tools repro repro-real corpus-repro acceptance uniref90-readiness clean

export PATH := $(CURDIR)/cache/bin:$(PATH)

sync:
	uv sync --extra dev

lint:
	uv run ruff check gasregnet scripts tests
	uv run mypy --strict gasregnet scripts/fetch_assets.py scripts/build_test_fixtures.py scripts/build_profiles.py scripts/run_corpus_discovery.py

test:
	uv run pytest -q --cov=gasregnet

assets:
	uv run gasregnet fetch-assets --manifest configs/assets.yaml --force

profiles:
	uv run gasregnet build-profiles --config configs --out-dir data/profiles --manifest-out data/profiles/profiles.yaml

seed-databases:
	uv run gasregnet build-seed-databases --config configs/headline.yaml --out data/profiles/diamond_seeds

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

discover-motifs embed-proteins foldseek-search fit-bayesian build-views emit-provenance:
	@echo "$@ is a placeholder target"

check-tools:
	uv run gasregnet check-tools --out tools_resolved.yaml

repro:
	uv run snakemake -s workflows/sqlite_mode.smk --cores 1

repro-real:
	uv run snakemake -s workflows/full_discovery.smk --cores 1

corpus-repro: datasets index-datasets
	uv run snakemake -s workflows/corpus_discovery.smk --cores 1

acceptance:
	./scripts/run_acceptance.sh

uniref90-readiness:
	uv run gasregnet check-tools --out /tmp/gasregnet-tools-resolved.yaml
	uv run gasregnet build-seed-databases --config configs/headline.yaml --out data/profiles/diamond_seeds
	uv run pytest tests/unit/test_config.py tests/unit/test_schemas.py tests/unit/datasets/test_refseq.py tests/unit/datasets/test_corpus_store.py tests/unit/annotation/test_ncbi_taxonomy.py tests/unit/search/test_hmmer.py
	uv run snakemake -s workflows/corpus_discovery.smk --dry-run --cores 4

clean:
	rm -rf results/* cache .snakemake .pytest_cache .ruff_cache .mypy_cache
	touch results/.gitkeep
