
## Plan: Real Sharding, Real DIAMOND, Real Acceptance

Six remaining gaps from the previous round are all *integration* gaps — every component exists in isolation but the corpus path doesn't actually use them. The plan is structured as five phases: **E** finishes the partitioned-store cutover, **F** lands real shard splitting, **G** wires DIAMOND seed-rescue per shard, **H** pins reproducibility (taxdump checksum + biology review), **I** stands up a real acceptance harness. Phases E→F→G are sequenced (each unblocks the next); H runs in parallel; I gates the readiness claim.

**TL;DR.** Today the corpus path has new modules but old call-sites: `extract_refseq_neighborhoods` still calls per-genome DuckDB; the Snakemake DAG splits *Python jobs* but not *catalogs*; seed-rescue is still the zip-by-index CPython loop; the taxdump and the scaling probe are unpinned. Five sequenced phases close all six gaps with a single end-state acceptance run.

Saved to /memories/session/plan.md.

---

### Phase E — Cut over anchor detection and neighborhood extraction to the partitioned store

**Goal.** refseq.py and anchors.py read exclusively from `data/corpus_store/{proteins,features}/` via Hive-partitioned `parquet_scan(...)` with one shared DuckDB connection. The per-genome `<dataset_name>.duckdb` files become unreferenced by the corpus path (kept as a back-compat option for the SQLite headline only).

**Steps**

1. Audit current call-sites. Use `grep_search` for `duckdb.connect` and `out_db` inside refseq.py and anchors.py. Confirm which functions still use the per-genome DuckDB path: at minimum `_anchor_feature_row`, `_read_cds_features`, `_feature_metadata`, `_all_protein_metadata`, plus `summarize_refseq_catalog` and `query_refseq_catalog`. Each must be migrated or explicitly retained behind a `legacy=True` flag.
2. Centralize the partitioned-store reader. New module gasregnet/datasets/corpus_reader.py exposing `open_corpus_store(store_root: Path) -> CorpusStoreHandle` returning a dataclass with a single open `duckdb.DuckDBPyConnection`, two registered views (`proteins_view`, `features_view`) over `parquet_scan(.../proteins/**/*.parquet, hive_partitioning=true)` and `.../features/**/*.parquet`, and helper methods `fetch_protein_metadata(accessions: list[str]) -> pl.DataFrame`, `fetch_cds_features(seqid: str, dataset_name: str) -> pl.DataFrame`, `fetch_anchor_feature(accession: str, dataset_name: str) -> dict`. The handle is a context manager (closes the connection on exit).
3. Migrate `_anchor_feature_row` and `_read_cds_features` in refseq.py to take a `CorpusStoreHandle` parameter rather than a `db: Path`. Eliminate per-anchor `duckdb.connect(...)` calls. Hot loop in `extract_refseq_neighborhoods` becomes a single iteration over `anchor_hits`, each iteration making 1–2 SQL calls against the *shared* connection, with predicate pushdown on the Hive partition columns (`phylum`, `dataset_name`).
4. Migrate `_feature_metadata` and `_all_protein_metadata` in anchors.py similarly. The seed-rescue loop's `_all_protein_metadata(db)` becomes `handle.fetch_proteins_for_dataset(dataset_name, phylum)` returning a Polars frame; this stays under DIAMOND replacement scope (Phase G) but uses the partitioned store now.
5. Update `index_refseq_corpus` and `summarize_refseq_corpus` to operate over `data/corpus_store/datasets.parquet` (the per-catalog provenance table from D.1). The legacy `<dataset_name>.duckdb` files become orphan; document a deprecation note in refseq.py module docstring and stop referencing them from any function on the corpus critical path.
6. Update the CLI surface: cli.py `extract_neighborhoods_command`, `detect_anchors_command`, `detect_anchors_profile_command`, `summarize_refseq_corpus_command`, `query_refseq_corpus_command` all accept a new `--store data/corpus_store/` flag; if absent, fall back to the legacy per-genome manifest path *with a deprecation warning emitted via structlog* so the V1 SQLite reproducibility leg keeps working unchanged.
7. Add a Pandera-level guard: when `provenance_source == "refseq"`, `LociSchema` requires the metadata column `corpus_store_root` (new optional column) to be non-null. This forces all corpus-path consumers downstream to flow through the partitioned store, not the legacy DuckDB.
8. Tests:
   - Unit: tests/unit/datasets/test_corpus_reader.py — verify single-connection reuse, partition pushdown (assert that querying for one `phylum` value scans only one partition file via DuckDB's `EXPLAIN`).
   - Integration: tests/integration/test_refseq_corpus_path_uses_store.py — run `extract_refseq_neighborhoods` end-to-end on a 5-catalog fixture, monkey-patch `duckdb.connect` to count invocations, assert `connect` is called ≤ 2 times total (one for the store, one for any explicit lookups), down from O(anchors).

**Relevant files**
- refseq.py → migrate `_anchor_feature_row`, `_read_cds_features`, `extract_refseq_neighborhoods`
- anchors.py → migrate `_feature_metadata`, `_all_protein_metadata`
- New gasregnet/datasets/corpus_reader.py
- cli.py → `--store` flag with legacy fallback
- schemas.py `LociSchema` → optional `corpus_store_root`
- New tests under datasets and integration

**Verification**
1. `grep -n "duckdb.connect" refseq.py gasregnet/search/anchors.py` returns zero hits inside any function reachable from `extract_refseq_neighborhoods` or `detect_anchors_profile` (legacy code paths gated behind `legacy=True` are acceptable).
2. `uv run pytest tests/integration/test_refseq_corpus_path_uses_store.py` passes, asserting connect-count ≤ 2.
3. Wallclock probe: re-run the 100-genome smoke from the previous readiness gate; per-anchor connect overhead drops from majority of wall-clock to < 1% (capture via `python -X importtime` or a structured timing log emitted by `corpus_reader`).
4. SQLite headline (sqlite_mode.smk) still passes byte-identical T1 output — the legacy SQLite path is untouched.

**Decisions**
- Keep the legacy per-genome DuckDB path behind a deprecation warning rather than deleting outright; this preserves the V1 reproducibility story and gives downstream users one release cycle to migrate.
- The `CorpusStoreHandle` is intentionally not a singleton — Snakemake fan-out (Phase F) will instantiate one per shard process, by design.

**Further considerations**
1. *Should `query_refseq_catalog` and `query_refseq_corpus` (the user-facing search commands) also migrate now?* **Recommendation: yes, in Phase E** — keeping two read paths invites drift. ~1 h additional work.

---

### Phase F — Real shard splitting in the Snakemake DAG

**Goal.** corpus_discovery.smk enumerates concrete shards from the partitioned store at DAG-construction time and dispatches one Snakemake job per shard. End state: `snakemake --cores 32 --dry-run` on a 1000-catalog corpus prints O(shards) parallel `detect_anchors_shard` jobs (not one monolithic shell).

**Steps**

1. Define the shard manifest. New module gasregnet/datasets/sharding.py: function `enumerate_shards(store_root: Path, strategy: ShardStrategy) -> list[Shard]` where `Shard = {shard_id: str, phylum: str, dataset_names: list[str], n_proteins_estimate: int}`. Two strategies: `by_phylum` (one shard per phylum value found in `datasets.parquet`) and `by_n_genomes` (greedy bin-packing of catalogs into shards of ≤ N genomes each). Persist the resolved shard manifest to `data/corpus_store/shards.parquet` so Snakemake input/output paths are stable across re-runs.
2. Generate the shard manifest via a new CLI command `gasregnet enumerate-shards --store data/corpus_store --strategy by_phylum --out data/corpus_store/shards.parquet`. This becomes a `rule enumerate_shards` in the Snakemake DAG that all per-shard rules depend on.
3. Rewrite corpus_discovery.smk using Snakemake `checkpoints`. The pattern:
   - `rule enumerate_shards` writes `shards.parquet` (a checkpoint, since shard count isn't known until indexing completes).
   - `checkpoint shards` reads `shards.parquet` and returns the list of `shard_id` wildcards.
   - `rule detect_anchors_shard` produces `intermediate/anchor_hits.{shard_id}.parquet` (one job per `shard_id`).
   - `rule extract_neighborhoods_shard` produces `intermediate/loci.{shard_id}.parquet`, `intermediate/genes.{shard_id}.parquet`.
   - `rule gather_anchors`, `rule gather_neighborhoods` concatenate via `pl.concat([pl.scan_parquet(p) for p in inputs]).collect()` and validate against the appropriate schema.
   - The downstream `annotate → assign_roles → score → enrich → archetypes → report` chain runs once on the gathered outputs.
4. Add per-shard CLI entry points that the rules invoke: `gasregnet detect-anchors-shard --store ... --shard-id <id> --out <path>` reads the shard's `phylum` and `dataset_names` filter from `shards.parquet`, materializes a Polars predicate over the partitioned store, and runs HMMER + DIAMOND (Phase G) on that subset only. Same for `extract-neighborhoods-shard`.
5. Add the SLURM profile workflows/profiles/slurm/config.yaml with per-rule resource hints — `detect_anchors_shard` rule gets `--cpus-per-task=8 --mem=16G --time=02:00:00`; `extract_neighborhoods_shard` gets `--cpus-per-task=4 --mem=8G`. Document `snakemake --profile workflows/profiles/slurm` invocation in workflows/profiles/README.md.
6. Add a local-execution profile workflows/profiles/local/config.yaml tuning thread counts for a single workstation; default to `cores: 8`.
7. Surface sharding configuration in corpus_discovery.yaml: `sharding: {strategy: by_phylum, n_genomes_per_shard: 1000, max_shards: null}`. The `max_shards` field exists for the scaling probe — it lets the acceptance harness cap shards to a known value.
8. Tests:
   - Unit: tests/unit/datasets/test_sharding.py — `by_phylum` produces one shard per distinct phylum; `by_n_genomes` greedy packing respects the cap and assigns every dataset exactly once; deterministic ordering for reproducibility.
   - Integration: tests/integration/test_corpus_dag_fanout.py — synthesize a 5-phylum × 4-catalog-per-phylum store; `snakemake --dry-run` should print 5 `detect_anchors_shard` jobs (with `by_phylum`) or 2 jobs (with `by_n_genomes` and `n_genomes_per_shard=10`).

**Relevant files**
- New gasregnet/datasets/sharding.py
- corpus_discovery.smk → checkpoint-based per-shard DAG
- New workflows/profiles/slurm/config.yaml, workflows/profiles/local/config.yaml, workflows/profiles/README.md
- corpus_discovery.yaml → sharding block
- cli.py → `enumerate-shards`, `detect-anchors-shard`, `extract-neighborhoods-shard` subcommands
- New tests as listed

**Verification**
1. `uv run gasregnet enumerate-shards --store data/corpus_store --strategy by_phylum --out /tmp/shards.parquet`; inspect with `duckdb -c "SELECT * FROM '/tmp/shards.parquet'"` — one row per phylum, deterministic shard_ids.
2. `uv run snakemake -s corpus_discovery.smk --cores 8 --dry-run` on a 5-phylum fixture prints exactly 5 `detect_anchors_shard` jobs.
3. `uv run snakemake -s corpus_discovery.smk --cores 8` on the same fixture produces `intermediate/anchor_hits.{shard_id}.parquet` files concurrently (verify via job timestamps in log).
4. `uv run pytest tests/integration/test_corpus_dag_fanout.py` passes.
5. The downstream gathered `intermediate/loci.parquet` is byte-identical (modulo row order) to a single-shard reference run on the same fixture — proves correctness of the gather step.

**Decisions**
- Snakemake `checkpoint` (not static expand) because shard count is data-dependent — phylum coverage is unknown until taxonomy resolution completes.
- `by_phylum` is the default strategy; `by_n_genomes` is for very-imbalanced corpora (e.g. UniRef90 dominated by Pseudomonadota).
- Per-shard outputs are named `{...}.{shard_id}.parquet` not `{shard_id}/{...}.parquet` because the gather step is a single concat, not a directory walk.

**Further considerations**
1. *Cluster vs single-node default?* **Recommendation: ship the local profile as default** in corpus_discovery.smk's `default-profile` directive; document SLURM as an opt-in in the README. UniRef90-scale runs will need the SLURM profile, but most users start single-node.

---

### Phase G — DIAMOND seed-rescue, per shard

**Goal.** Replace `_seed_rescue_hits` (CPython zip-by-index against every protein in a catalog) with a per-shard DIAMOND `--ultra-sensitive` invocation. End state: seed-rescue wallclock on a 1000-catalog shard drops from "minutes per catalog" to "seconds per shard total".

**Steps**

1. Build seed DIAMOND databases at indexing time, not at search time. New CLI subcommand `gasregnet build-seed-databases --config headline.yaml --out data/profiles/diamond_seeds/`. For each `(analyte, anchor_family)` pair, runs `diamond makedb --in <seeds.faa> --db data/profiles/diamond_seeds/<analyte>__<family>.dmnd`. Persist version + DB hash in profiles.yaml under each profile row's new `seed_diamond_db` key.
2. New module gasregnet/search/seed_rescue.py with public function `seed_rescue_for_shard(*, shard: Shard, store: CorpusStoreHandle, seeds_diamond_dir: Path, analytes: list[AnalyteConfig], identity_threshold: float = 0.30, coverage_threshold: float = 0.50, e_value: float = 1e-10, threads: int = 8) -> pl.DataFrame`. Implementation:
   - Materialize the shard's proteins to a temp FASTA via `store.fetch_proteins_for_shard(shard)` writing through `gasregnet.io.fasta.write_fasta`.
   - For each `(analyte, anchor_family)` pair, run `diamond blastp` with the per-pair DB as subject (small) and the shard FASTA as query (large), using `--ultra-sensitive --outfmt 6 ...`. This is the DIAMOND-recommended direction (many queries, small DB).
   - Parse via existing `parse_diamond_output`, filter by `(percent_identity, qcovhsp)` thresholds, keep best hit per (query, anchor_family).
   - Apply the family guard from `_passes_family_guard` as a vectorized Polars filter against `product`/`gene`/`locus_tag` (already audited, was lambda-based).
   - Return a frame conformant to `AnchorHitsSchema` with `evidence_type = "seed_back_confirmed"`.
3. Replace the call in anchors.py `detect_anchors_profile`: when `--store` is set, call `seed_rescue_for_shard` per shard rather than `_seed_rescue_hits` per catalog. Preserve `_seed_rescue_hits` behind the `legacy=True` flag for the SQLite path's regression tests; mark deprecated.
4. Wire into the per-shard rule. `gasregnet detect-anchors-shard` (from Phase F step 4) becomes the union of HMMER profile hits + DIAMOND seed-rescue hits for that shard, deduplicated on `(dataset_name, analyte, anchor_family, protein_accession, locus_tag)`.
5. Vectorize the family guard. In anchors.py, replace the `pl.struct(...).map_elements(_passes_family_guard, ...)` lambda with `pl.col("product").str.to_lowercase().str.contains_any(terms) | pl.col("gene")...`. This kills the audit's Pf8 finding inside the same phase.
6. Stream HMMER inputs. In hmmer.py, replace `list(sequence_file)` with the streaming `SequenceFile` iterator passed to `pyhmmer.hmmer.hmmsearch`. At UniRef90 50 GB FASTA scale this is the difference between OOM and not (audit Pf7).
7. Tests:
   - Unit: tests/unit/search/test_seed_rescue.py — fixture with one CoxL homolog and one decoy, assert DIAMOND call structure (mock subprocess), assert AnchorHitsSchema conformance.
   - Integration: tests/integration/test_seed_rescue_recovers_partial_homologs.py — real DIAMOND run on a small fixture protein set including a 60%-identity CoxL distant homolog that HMMER profile misses; assert seed-rescue recovers it.
   - Performance gate: tests/integration/test_seed_rescue_scales.py — synthetic shard of 10K proteins; assert wallclock < 30 s on 4 cores. Skipped if `diamond` binary unavailable.
8. Ensure DIAMOND is in the `check-tools` preflight (Phase H step 5 below); fail loudly with actionable error when missing.

**Relevant files**
- New gasregnet/search/seed_rescue.py
- anchors.py → swap `_seed_rescue_hits` for `seed_rescue_for_shard`; vectorize family guard
- hmmer.py → stream SequenceFile
- cli.py → `build-seed-databases` subcommand
- profiles.yaml → `seed_diamond_db` column per profile row
- New tests as listed

**Verification**
1. `uv run gasregnet build-seed-databases --config headline.yaml --out data/profiles/diamond_seeds/` produces one `.dmnd` per `(analyte, family)` pair; `ls data/profiles/diamond_seeds/` shows expected count.
2. `uv run pytest tests/integration/test_seed_rescue_recovers_partial_homologs.py` passes.
3. `uv run pytest tests/integration/test_seed_rescue_scales.py` passes (gated on local DIAMOND availability).
4. Wallclock on the 100-genome readiness smoke: total seed-rescue time ≤ 5% of `detect_anchors_shard` wallclock (was the dominant cost pre-fix). Capture via structured timing log in `seed_rescue_for_shard`.
5. Recall regression: on the V2 benchmark, seed-rescue+HMMER union recovers ≥ as many anchors as the old `_seed_rescue_hits`+HMMER union (no recall loss from migration).

**Decisions**
- DIAMOND as subject-DB-small / query-FASTA-large is the standard "search many sequences against a small library" pattern; reverses the direction the legacy code used.
- Per-pair `.dmnd` files (not one monolithic DB) because anchor families are independent and we want per-family score thresholds.
- Default `--ultra-sensitive` not `--very-sensitive` because seed-rescue is exactly the partial-homolog recall scenario where ultra-sensitive's MMseqs-comparable recall justifies its 2-3× wallclock cost.

**Further considerations**
1. *MMseqs2 as a future ablation?* Out of scope here, but the `seed_rescue_for_shard` interface accepts a `tool: Literal["diamond","mmseqs2"] = "diamond"` parameter as a forward-compatible hook.

---

### Phase H — Reproducibility pinning + biology review

**Goal.** Every external input the corpus path consumes is checksum-pinned in the manifest; CO and CN curated content is reviewed by a domain expert before manuscript claims. End state: two consecutive runs of the readiness gate produce identical `manifest.json` `run_hash`; biology review sign-off is filed alongside seeds and benchmarks.

**Steps**

1. Pin the NCBI taxdump asset. Update assets.yaml entry for `ncbi_taxdump`: add `sha256: <pinned hash>`, `version: <YYYY-MM-DD release tag>`, `source_url: ...`. Compute the hash from the actually-shipped taxdump (`shasum -a 256 data/external/ncbi_taxonomy/taxdump.tar.gz`); commit the value. Update assets.py `Downloader` to verify the SHA-256 after download and *fail* if it doesn't match (not warn). Same for the taxdump expansion: add `sha256` for `nodes.dmp` and `names.dmp` after un-tarring.
2. Pin DIAMOND seed databases. Each `.dmnd` is non-deterministic across DIAMOND versions; pin via the *DIAMOND version + seed FASTA hash* tuple, not the `.dmnd` byte hash. Add `seed_diamond_db_version: <diamond version>` and `seed_fasta_sha256: <hash>` columns to profiles.yaml per profile row. Re-run `build-seed-databases` whenever either changes.
3. Implement `gasregnet check-tools`. New CLI subcommand that probes each external binary (`hmmsearch -h`, `diamond version`, optionally `foldseek version`, `mafft --version`, `meme -version`) and writes `tools_resolved.yaml` to the run directory with `{binary: {version, path, sha256_of_binary}}`. Wire into manifest.py: every CLI command persists `external_versions` from `tools_resolved.yaml` into the manifest. Closes audit T5.
4. Make `manifest.json` `run_hash` stable across runs. Audit `build_manifest` to ensure every input file referenced (config, seed FASTA, profile HMM, taxdump, benchmark CSV) is hashed; add the missing ones (audit T4). Add a regression test tests/integration/test_manifest_run_hash_stable.py that runs the SQLite headline twice and asserts identical `run_hash`.
5. Biology review handoff. Create `data/seeds/REVIEW_LOG.md` and `data/benchmarks/REVIEW_LOG.md` (these are review records, not docs of code changes — explicitly within the user's "do not document changes" exemption since they're scientific-process artifacts). Each entry: reviewer name+affiliation, review date, files reviewed, scope (`co_coxL_seeds.faa`, `co_cooS_seeds.faa`, `regulators_v2_cn_curated.csv`, `known_cn_organisms.csv`), changes requested, signed-off git SHA. Block merging the V2 manuscript draft until both review logs have at least one signed entry from a non-author domain expert. Pre-fill a review-prompt template that asks specifically: "(1) does the CoxL seed set exclude CooS contamination at the family level; (2) does the CooS seed set exclude bacterial vs archaeal CdhA confounding; (3) is the HCN biosynthesis benchmark complete enough to support the V2 chemistry-partition claim".
6. Pre-flight checksum probe. Add a Snakemake `rule preflight` that runs at DAG entry: verifies all pinned assets match their declared SHA-256 and that all required binaries pass `check-tools`. Fails the workflow with an actionable error before any expensive work starts.

**Relevant files**
- assets.yaml → pinned `sha256`, `version`, `source_url` for taxdump
- assets.py → enforce SHA-256 verification on download
- profiles.yaml → `seed_diamond_db_version`, `seed_fasta_sha256` columns
- cli.py → `check-tools` subcommand
- manifest.py → consume `tools_resolved.yaml`; expand input hashing
- New `data/seeds/REVIEW_LOG.md`, `data/benchmarks/REVIEW_LOG.md` (scientific-process artifacts)
- corpus_discovery.smk → `rule preflight` at DAG entry
- New test tests/integration/test_manifest_run_hash_stable.py

**Verification**
1. `uv run gasregnet fetch-assets --manifest assets.yaml --force` — modify the local taxdump byte and re-run; the second invocation fails with a SHA-256 mismatch error pointing at the specific asset.
2. `uv run gasregnet check-tools --out /tmp/tools_resolved.yaml`; inspect — every required binary has version + path + SHA-256.
3. `uv run pytest tests/integration/test_manifest_run_hash_stable.py` passes (two consecutive runs → identical `run_hash`).
4. `data/seeds/REVIEW_LOG.md` has at least one signed entry referring to a specific git SHA before any UniRef90-scale claim is made.
5. `uv run snakemake -s corpus_discovery.smk --cores 1 -n` prints `preflight` as the first job; manually delete a pinned binary and confirm the workflow refuses to start.

**Decisions**
- SHA-256 *failure*, not warning, for asset verification — silent reproducibility drift is the entire problem we're trying to prevent.
- Biology review gate is a process gate, not an automated test; relying on a CODEOWNERS-style approval for seeds and benchmarks PRs is the lighter-weight alternative if a formal review log is overkill — call this out in the followup.
- DIAMOND DBs are pinned by `(version, seed_fasta_sha256)`, not by `.dmnd` byte hash, because DIAMOND DB format is non-deterministic across versions.

**Further considerations**
1. *Where does the biology reviewer come from?* **Recommendation: explicit named expert** rather than a community PR review. CO/anaerobic-microbiology and HCN-biosynthesis-microbiology are specialized; one named reviewer per analyte. Identify before Phase B's curation work begins (already on the audit critical path).

---

### Phase I — Real acceptance harness (10K-genome scaling probe)

**Goal.** Replace the readiness-gate scaffolding with a real, runnable, end-to-end probe that exercises 10,000 genomes through the full corpus path and asserts wallclock + memory + correctness budgets. Pass = "we believe the code is ready to launch on UniRef90". Fail = "stop and fix before launch".

**Steps**

1. Generate the 10K synthetic corpus. New script scripts/generate_scaling_corpus.py: produces 10K synthetic NCBI-style assembly catalogs into `data/scaling_probe/<dataset_name>/` (each containing `protein.faa.gz`, `genomic.gff.gz`, `assembly_report.txt` with a planted `# Taxid:` line). Genome size: 1000 proteins each (small but realistic). Phylum distribution: roughly mirror UniRef90 — Pseudomonadota 35%, Bacillota 18%, Bacteroidota 12%, Actinomycetota 10%, ~25 other phyla for the remainder. Plant a known number of CoxL, CooS, HcnA, FNR proteins (say 50 of each across the corpus) at known positions for recall validation.
2. Generate a synthetic taxdump that resolves the planted Taxids. Script: scripts/generate_scaling_taxdump.py emits a minimal `nodes.dmp`/`names.dmp` covering exactly the 10K corpus's Taxids. Avoids needing to ship the real taxdump in CI.
3. Define the acceptance budgets. New file tests/integration/scaling_probe/budgets.yaml with version-controlled targets:
   - `wallclock_seconds: 3600` (1 h on 8 cores; arithmetic: 10K genomes × ~10 anchors × ~10 ms per anchor SQL = 1000 s ideal + DIAMOND seed-rescue ~30 s/shard × ~25 phylum shards = 750 s, plus annotate/score/enrich amortized; 1 h is generous headroom)
   - `peak_memory_gb: 32`
   - `connect_count_max: 100` (one per shard process, with margin)
   - `recall_floor: {coxL: 0.95, cooS: 0.95, hcnA: 0.90, fnr: 0.90}` — 95% of planted CoxL and CooS proteins must appear in `anchor_hits.parquet`
   - `partition_distinct_chemistries_min: 4` (heme, iron_sulfur_4Fe4S, flavin, none — minimum diversity for Figure 4 to be non-degenerate)
   - `taxon_id_zero_count_max: 0` — schema regression check
4. Build the acceptance test tests/integration/scaling_probe/test_acceptance.py. Marked `@pytest.mark.scaling` so it's opt-in. Steps:
   - Build the synthetic corpus (cached in `data/scaling_probe/` with checksum; regenerate if missing).
   - Run `gasregnet index-refseq-corpus --store data/scaling_probe/store/ --taxonomy-db data/scaling_probe/taxonomy.duckdb`.
   - Run `snakemake -s corpus_discovery.smk --cores 8 --profile workflows/profiles/local`.
   - Capture wallclock from `time` + memory from `/usr/bin/time -v` (or psutil polling).
   - Capture connect-count via a structlog filter.
   - Read the final `loci.parquet`, `anchor_hits.parquet`, `enrichment.parquet`.
   - Assert each budget. Emit a structured report `tests/integration/scaling_probe/last_run_report.json` with the actuals.
5. CI integration. The acceptance test does not run in default `pytest`; it runs in a dedicated GitHub Actions job (or Make target `make acceptance`) on a self-hosted runner with adequate resources. Add `scripts/run_acceptance.sh` that wraps the invocation and uploads the report.
6. Failure-mode dashboard. The acceptance test on failure must produce actionable output: which budget failed, by how much, with a pointer to the offending phase (e.g. "wallclock exceeded budget by 1.4×; per-rule timing in `intermediate/snakemake_timings.tsv` shows `detect_anchors_shard` is the bottleneck"). This means structured per-rule timing logs from Snakemake (`benchmark` directive on each per-shard rule) and structured per-function timing inside `corpus_reader`.
7. Pre-launch UniRef90 dry-run script. New scripts/uniref90_dryrun.sh that takes the *real* UniRef90 release path, runs `enumerate-shards --strategy by_phylum --max-shards 5`, and executes only the first 5 shards through `detect-anchors-shard`. Catches issues that surface only on real heterogeneous data before committing to a multi-day full-corpus run. Document the workflow in workflows/profiles/README.md.

**Relevant files**
- New scripts/generate_scaling_corpus.py
- New scripts/generate_scaling_taxdump.py
- New tests/integration/scaling_probe/budgets.yaml
- New tests/integration/scaling_probe/test_acceptance.py
- New scripts/run_acceptance.sh
- New scripts/uniref90_dryrun.sh
- Makefile → new target `acceptance:` invoking `scripts/run_acceptance.sh`
- corpus_discovery.smk → `benchmark:` directives on per-shard rules

**Verification**
1. `make acceptance` on a workstation with 8 cores, 32 GB RAM completes within 1 hour, peak RSS < 32 GB; produces `tests/integration/scaling_probe/last_run_report.json` with all budgets in `pass` state.
2. Recall on planted proteins: ≥ 95% for CoxL, CooS; ≥ 90% for HcnA, FNR. Re-injecting noise (random sequence shuffling at 30% identity) should drop recall below the floor — confirms the test isn't trivially passing.
3. Connect-count budget is met (≤ 100 across the entire 10K-genome run) — proves Phase E migration held under fan-out.
4. Re-running the acceptance test produces identical `run_hash` in the final manifest — proves Phase H pinning held.
5. Deliberate-failure injection: revert one Phase E migration commit; acceptance test fails with `connect_count_max` exceeded and clearly identifies `extract_refseq_neighborhoods` as the offender.

**Decisions**
- 10K genomes (not 100, not 100K): 100 is too small to expose scale defects; 100K requires a CI runner not many people have. 10K fits in 1 h × 8 cores and exposes every defect class.
- Synthetic corpus, not a real RefSeq subset: real subsets bring real download/network/checksum surface area into CI; synthetic with planted positives gives deterministic recall floors.
- Budgets are version-controlled YAML so reviewers can see exactly what "ready" means and so tightening budgets is a deliberate code change.
- Pre-launch UniRef90 dry-run is a separate, manual step — not in CI, but documented for the operator who actually launches the production run.

**Further considerations**
1. *Should the acceptance test run on PRs or only on main?* **Recommendation: nightly on main + manual `workflow_dispatch` trigger**. Per-PR is too expensive at 1 h. Block merges to `release/*` branches on a green nightly.
2. *What's the failure budget?* I.e. how many flaky-CI retries before declaring not-ready. **Recommendation: zero retries on the wallclock budget** (deterministic failure), up to 3 retries on the recall floor (DIAMOND scoring has minor stochastic variation only with `--ultra-sensitive` heuristic seed sampling).

---

### Sequencing & ownership

```
Phase E ─────────────► Phase F ─────────────► Phase G ─────────────► Phase I
                                                                       ▲
Phase H (parallel) ────────────────────────────────────────────────────┘
```

- **Phase E** unblocks F (sharding needs the partitioned-store reader) and G (DIAMOND seed-rescue per shard needs the store handle).
- **Phase F** unblocks G (DIAMOND runs *per shard*, not per catalog).
- **Phase G** unblocks I (acceptance probe exercises seed-rescue at scale).
- **Phase H** is parallel-safe with E/F/G. The biology review process should kick off *now* in parallel with engineering work.
- **Phase I** is the readiness gate; it must pass before any UniRef90 launch.

Suggested agent allocation:
- Agent 1: Phase E → Phase G (storage + tool integration)
- Agent 2: Phase F (workflow / Snakemake) — starts after E lands the `CorpusStoreHandle` interface
- Agent 3: Phase H (reproducibility) — fully parallel
- Agent 4 + named domain reviewer: Phase H step 5 (biology review)
- Agent 5: Phase I — starts after E, F, G all land

### What I deliberately excluded

- ESM/Foldseek/PyMC integrations from the audit — orthogonal to scale readiness.
- MMseqs2 as alternative seed-rescue — `seed_rescue_for_shard` has the hook, implementation is post-launch.
- GTDB lineage as alternative to NCBI Taxonomy — same hook (`provenance_source`-conditional), post-launch.
- Foldseek-based `structural_plausibility_score` wiring — separate work-stream.

If you approve, the same handoff pattern as last round: Agents 1, 3, 4 start immediately on E, H, and biology-review-prep; Agent 2 picks up F as soon as E lands the `CorpusStoreHandle` contract; Agent 5 builds the acceptance harness in parallel with G so the moment G lands we have a real readiness signal. Want me to persist this plan and proceed, or revise any phase first?Continued with Autopilot: Task not yet completePersisting the plan to session memory now, then completing. 

Read memory [](file:///Users/filiprumenovski/Library/Application%20Support/Code/User/workspaceStorage/6797a959b0679f66e69e1e5982d32470/GitHub.copilot-chat/memory-tool/memories/OTFhNmY5YzMtMmI2NS00ZTJkLWJiNTAtMTkzYmQ4MGUzZjcz/plan.md)

Created memory file [](file:///Users/filiprumenovski/Library/Application%20Support/Code/User/workspaceStorage/6797a959b0679f66e69e1e5982d32470/GitHub.copilot-chat/memory-tool/memories/OTFhNmY5YzMtMmI2NS00ZTJkLWJiNTAtMTkzYmQ4MGUzZjcz/plan.md)

Plan persisted to /memories/session/plan.md. Five phases (E storage cutover → F real shard splitting → G DIAMOND seed-rescue → I acceptance harness, with H reproducibility pinning + biology review running in parallel) close the six remaining UniRef90-readiness gaps. Full plan was presented in the previous turn. 

