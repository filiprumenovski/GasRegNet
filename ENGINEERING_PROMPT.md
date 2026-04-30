# GasRegNet Engineering Prompt

Ground-truth specification for coding agents. Pair with `ARCHITECTURE.md`
(scientific scope) and `PREPRINT_PLAN.md` (publication target). Read all three
before writing code.

## 0. How to use this document

Every agent reads sections 1 through 7 in full before writing code. Sections
8 onward are unit specifications. Pick exactly one task at a time from
section 11. Implement against the data contracts in section 6 and the
conventions in section 5. Run the tests in section 12 before declaring the
task done. Boundary contracts between modules are enforced by Pandera schemas,
not by trust.

This codebase will be developed by multiple parallel agents. Module boundaries
are designed so that two agents working on different modules at the same time
should never produce conflicting changes to the same file. If you find
yourself wanting to modify a file outside your assigned module, stop and
record the cross-cut as an issue rather than editing.

If anything in this document conflicts with your own judgment about the
"right" abstraction, defer to this document. The user's domain expert wrote
it. Your job is to build, not to redesign.

## 1. Mission

Build a Python package, `gasregnet`, that takes bacterial genome and protein
data as input and produces ranked, statistically validated candidate gas-sensing
transcriptional regulators as output. The package must be runnable in three
modes: a fast SQLite mode that consumes EFI-GNT exports, a curated-panel mode
that consumes user-supplied FASTA and GFF, and a full-discovery mode that runs
DIAMOND searches over a large proteome database. All three modes share the
same downstream scoring, archetype, enrichment, and reporting modules.

The output of the package is not raw hits. It is a publication-ready set of
CSV tables, Parquet intermediates, and figures supporting three falsifiable
claims: that the method recovers known canonical bacterial gas-responsive
regulators on a versioned benchmark; that it nominates novel candidates that
pass matched-control enrichment at FDR < 0.05; and that CO and HCN regulators
partition into chemically distinct architectures. The package must produce
all assets needed to draft a bioRxiv preprint. Reproducibility is a
first-class requirement: every run is deterministic given pinned configs,
and a `make repro` target reconstructs all main figures from a clean clone.

The user is a single researcher orchestrating multiple LLM coding agents in
parallel. Code must be readable, modular, tested, and conservative about new
dependencies. Treat the existing scripts in `GNN_analysis_tool/` as historical
context, not as a base to extend.

## 2. Non-negotiables

These are correctness requirements, not preferences.

1. Every public function returns a typed object validated by a Pandera schema
   or a Pydantic model. No untyped DataFrames cross module boundaries.
2. Every CLI command accepts `--config` pointing to a YAML file or directory
   and `--out` pointing to a results directory. Configs are loaded, validated
   against a Pydantic schema, and echoed to `<out>/config.resolved.yaml` at
   the start of every run.
3. Every results directory contains `manifest.json` with run hash, package
   version, config hashes, input data hashes, external binary versions, seed,
   wall clock, and timestamp.
4. Polars is the default DataFrame library. Pandas is allowed only at
   integration boundaries with libraries that require it.
5. Snakemake orchestrates multi-step pipelines. Do not introduce ad-hoc shell
   scripts that bypass Snakemake.
6. Logging uses `structlog` at `INFO` by default and `DEBUG` under
   `--verbose`. No `print()` calls outside of CLI user-facing output.
7. All randomness is seeded. The default seed is read from config. Default
   value is `20260429`.
8. No internet access at runtime in the SQLite or curated-panel modes.
   Full-discovery mode may use external APIs only when explicitly enabled in
   config.
9. Tests live next to the code they exercise. Coverage on `gasregnet/` is
   enforced at 80% minimum for any module marked P1 or P2 in the task graph.
10. No silent fallbacks. If a config field is missing, raise. If a benchmark
    protein is not in the corpus, log a warning, continue, and record the
    miss in the manifest.

## 3. Tech stack and pinning

Pin exact versions in `pyproject.toml`. Lock with `uv` and commit `uv.lock`.

```text
Python:               3.11.9
Package manager:      uv 0.5.x
Build system:         hatchling
DataFrame:            polars 1.x
Schema validation:    pandera[polars] 0.23.x
Config validation:    pydantic 2.x
Workflow:             snakemake 8.x
Storage:              duckdb 1.x, pyarrow 17.x
Logging:              structlog 24.x
CLI:                  typer 0.12.x
Statistics:           scipy 1.13.x, statsmodels 0.14.x
Plotting:             matplotlib 3.9.x, seaborn 0.13.x
Bio:                  biopython 1.84
HMM:                  pyhmmer 0.10.x
Search:               diamond binary 2.1.10 (external)
Cluster:              mmseqs2 binary 15-6f452 (external)
Domains:              interproscan 5.x (external, optional)
Structure:            none in-package; ingest AlphaFold-Server outputs only
Linting:              ruff 0.7.x
Typing:               mypy 1.11.x
Testing:              pytest 8.x, hypothesis 6.x, pytest-cov 5.x
```

External binaries are not pinned in Python. Their versions are recorded in
`manifest.json` at runtime by parsing `--version`.

## 4. Repository layout

The exact tree. Create empty `__init__.py` where appropriate.

```text
GasRegNet/
  ARCHITECTURE.md                    # do not edit
  PREPRINT_PLAN.md                   # do not edit
  ENGINEERING_PROMPT.md              # this file; do not edit
  CLAIMS.md                          # human-authored in T1
  README.md                          # produced in T2
  CITATION.cff                       # produced at v0.1.0 tag
  pyproject.toml
  uv.lock
  Makefile
  configs/
    analytes/
      co.yaml
      cn.yaml
    regulator_families.yaml
    sensory_domains.yaml
    scoring.yaml
    benchmarks.yaml
    headline.yaml
  data/
    seeds/
      co_anchor_seeds.faa
      cn_anchor_seeds.faa
    benchmarks/
      benchmark_v1.csv               # produced in T0
    references/
      known_co_organisms.csv
      known_cn_organisms.csv
    controls/                        # populated at runtime
  gasregnet/
    __init__.py
    cli.py
    config.py
    schemas.py
    errors.py
    paths.py
    logging.py
    hashing.py
    manifest.py
    io/
      sqlite_efi.py
      fasta.py
      gff.py
      parquet.py
    search/
      diamond.py
      hmmer.py
      mmseqs.py
    annotation/
      domains.py
      regulators.py
      taxonomy.py
      ecology.py
    neighborhoods/
      retrieve.py
      encode.py
      operons.py
    scoring/
      loci.py
      candidates.py
      enrichment.py
      controls.py
    archetypes/
      cluster.py
      diagrams.py
    structure/
      msa.py
      alphafold.py
      pdb.py
    reports/
      figures.py
      tables.py
      captions.py
      gasregnet.mplstyle
  workflows/
    sqlite_mode.smk
    diamond_mode.smk
    full_discovery.smk
  tests/
    unit/
      io/
      search/
      annotation/
      neighborhoods/
      scoring/
      archetypes/
      structure/
      reports/
    integration/
      test_sqlite_end_to_end.py
      test_benchmark_recovery.py
      test_partition_claim.py
    fixtures/
      mini_efi.sqlite
      mini_co_seeds.faa
      mini_cn_seeds.faa
      mini_genomes.tsv
      mini_alphafold/
  results/                           # gitignored; runtime outputs
    .gitkeep
  scripts/
    bootstrap_dev.sh
    fetch_external_binaries.sh
    build_test_fixtures.py
  .github/
    workflows/
      ci.yml
```

## 5. Conventions

Type every function. `from __future__ import annotations` at the top of every
module. Run `mypy --strict` on `gasregnet/` and resolve all errors before
declaring a task done.

Logging is structured. Every log line carries `module`, `function`, and any
IDs in scope.

```python
import structlog
log = structlog.get_logger(__name__)
log.info("loaded loci", n_loci=42, analyte="CO")
```

Errors are raised as typed exceptions defined in `gasregnet/errors.py`.
Define `ConfigError`, `SchemaError`, `MissingInputError`, `BenchmarkMissError`,
`ScoringError`, `EnrichmentError`. Do not raise bare `Exception`. Do not
silently swallow exceptions; catch typed, log with context, re-raise or
convert.

Paths are always `pathlib.Path`. Convert to `str` only at the boundary with
libraries that require it.

Configs flow as Pydantic models. The CLI loads YAML, validates, and passes
the model. Do not pass dicts.

Determinism. Every function that uses randomness accepts `seed: int` with no
default. The CLI supplies the seed from config.

Output. Every run writes to a single results directory containing
`manifest.json`, `config.resolved.yaml`, `logs/run.jsonl`, and subdirectories
`tables/`, `figures/`, `intermediate/`. Never write outside `--out`.

Naming. Snake case for modules and functions. Pascal case for classes. Lower
snake case for table columns. Use bioinformatics canon abbreviations
(`pfam_id`, `gff`, `faa`, `msa`, `pdb`, `fdr`, `roc`, `pr`). Spell out
otherwise.

Imports. Standard library, then third party, then `gasregnet.*`. Sort within
each block. Ruff enforces this.

## 6. Data contracts

Every cross-module table is defined as a Pandera schema in
`gasregnet/schemas.py`. Schemas are the single source of truth. Tables on
disk are Parquet, validated on read and on write.

Field types use Polars dtypes.

### `LociSchema`

```text
locus_id:                        Utf8       primary key, format f"{analyte}_{cluster_id}_{anchor_accession}"
analyte:                         Utf8       in {"CO", "CN"}
anchor_accession:                Utf8       UniProt or NCBI protein accession
anchor_family:                   Utf8       e.g., "coxL", "cydA"
organism:                        Utf8
taxon_id:                        Int64      NCBI taxon
cluster_id:                      Int32
contig_id:                       Utf8
window_size:                     Int32
is_boundary_truncated:           Boolean
marker_genes_present:            List[Utf8]
accessory_genes_present:         List[Utf8]
locus_score:                     Float64
locus_confidence:                Utf8       in {"high", "medium", "low", "control"}
taxonomic_context_score:         Float64
operon_integrity_score:          Float64
created_at:                      Datetime(time_unit="us")
```

### `GenesSchema`

```text
locus_id:                        Utf8       FK to LociSchema
gene_accession:                  Utf8
relative_index:                  Int32      0 = anchor, negative upstream
relative_start:                  Int64
relative_stop:                   Int64
strand:                          Utf8       in {"+", "-"}
product_description:             Utf8
pfam_ids:                        List[Utf8]
pfam_descriptions:               List[Utf8]
interpro_ids:                    List[Utf8]
interpro_descriptions:           List[Utf8]
functional_class:                Utf8       in {"anchor", "regulator", "metabolic", "transporter", "unknown"}
regulator_class:                 Utf8       in {"one_component", "two_component_rr", "two_component_hk", "sigma", "anti_sigma", "antiterminator", "none"}
sensory_domains:                 List[Utf8]
is_anchor:                       Boolean
is_regulator_candidate:          Boolean
```

Cross-field rule: `relative_index == 0` implies `is_anchor == True`.

### `RegulatorCandidatesSchema`

```text
candidate_id:                    Utf8       primary key
analyte:                         Utf8
locus_id:                        Utf8       FK to LociSchema
gene_accession:                  Utf8
organism:                        Utf8
cluster_id:                      Int32
relative_index:                  Int32
distance_nt:                     Int64
position:                        Utf8       in {"upstream", "downstream", "internal"}
strand:                          Utf8
regulator_class:                 Utf8
dna_binding_domains:             List[Utf8]
sensory_domains:                 List[Utf8]
pfam_ids:                        List[Utf8]
interpro_ids:                    List[Utf8]
archetype_id:                    Utf8       FK to ArchetypesSchema, nullable
locus_score:                     Float64
regulator_domain_score:          Float64
sensory_domain_score:            Float64
proximity_score:                 Float64
archetype_conservation_score:    Float64
enrichment_score:                Float64
taxonomic_breadth_score:         Float64
structural_plausibility_score:   Float64    nullable
candidate_score:                 Float64
candidate_score_q:               Float64    nullable, in [0, 1]
rationale:                       Utf8
```

### `ArchetypesSchema`

```text
archetype_id:                    Utf8       primary key
analyte:                         Utf8
cluster_id:                      Int32
architecture_string:             Utf8       e.g., "[-3:HK:PAS][-2:RR:REC-HTH][-1:cydB][0:cydA][+1:cydX]"
n_loci:                          Int32
n_taxa:                          Int32
representative_locus_id:         Utf8       FK to LociSchema
dominant_regulator_class:        Utf8
dominant_anchor_structure:       Utf8
mean_locus_score:                Float64
mean_candidate_score:            Float64
```

### `EnrichmentResultsSchema`

```text
analyte:                         Utf8
feature_type:                    Utf8       in {"regulator_family", "sensory_domain", "regulator_class", "archetype"}
feature_name:                    Utf8
case_definition:                 Utf8
control_definition:              Utf8
n_case_with_feature:             Int64
n_case_without_feature:          Int64
n_control_with_feature:          Int64
n_control_without_feature:       Int64
odds_ratio:                      Float64
p_value:                         Float64
q_value:                         Float64
interpretation:                  Utf8
```

### `BenchmarkSchema`

```text
benchmark_id:                    Utf8       primary key
analyte:                         Utf8       in {"CO", "CN", "negative_control"}
protein_name:                    Utf8       e.g., "CooA", "RcoM-1", "MpaR"
uniprot_accession:               Utf8
organism:                        Utf8
taxon_id:                        Int64
anchor_family:                   Utf8       nullable for negative controls
expected_regulator_class:        Utf8
expected_sensory_domains:        List[Utf8]
sensing_evidence_class:          Utf8       in {"direct", "indirect", "inferred", "none"}
pmid:                            List[Utf8]
notes:                           Utf8
```

## 7. Configuration files

All configs validate against Pydantic models in `gasregnet/config.py`. Below
are the canonical instances.

### `configs/analytes/co.yaml`

```yaml
analyte: CO
display_name: "carbon monoxide"
anchor_seeds: data/seeds/co_anchor_seeds.faa
anchor_families:
  - name: coxL
    pfam_required: ["PF02738"]
    pfam_supporting: ["PF03450", "PF01799"]
    role: primary
  - name: coxM
    pfam_required: ["PF03450"]
    role: accessory
  - name: coxS
    pfam_required: ["PF01799"]
    role: accessory
window_genes: 10
known_organisms_table: data/references/known_co_organisms.csv
expected_sensory_chemistry:
  - heme
  - iron_sulfur
seed: 20260429
```

### `configs/analytes/cn.yaml`

```yaml
analyte: CN
display_name: "hydrogen cyanide"
anchor_seeds: data/seeds/cn_anchor_seeds.faa
anchor_families:
  - name: cydA
    pfam_required: ["PF01654"]
    role: primary
  - name: cydB
    pfam_required: ["PF02322"]
    role: primary
  - name: cydX
    pfam_required: ["PF13456"]
    role: accessory
window_genes: 10
known_organisms_table: data/references/known_cn_organisms.csv
expected_sensory_chemistry:
  - cysteine_rich
  - flavin
  - heme
seed: 20260429
```

### `configs/regulator_families.yaml`

Maps Pfam and InterPro IDs to regulator families and classes. Required entries
include FNR/CRP (PF00027 cNMP binding plus HTH PF00325), LysR (PF00126 plus
PF03466), GntR (PF00392), AraC (PF12833), MarR (PF01047), TetR (PF00440),
MerR (PF00376), LacI (PF00356), ArsR/SmtB (PF01022), LuxR/GerE (PF00196),
XRE (PF01381), two-component receiver (PF00072), histidine kinase (PF00512),
HATPase (PF02518), sigma-70 (PF04542, PF04545). Each entry has `family`,
`class`, `pfam_required`, `pfam_optional`.

### `configs/sensory_domains.yaml`

Maps Pfam and InterPro IDs to sensory chemistries. PAS (PF00989, PF13426,
PF14598), GAF (PF01590), Cache (PF02743, PF13185), HAMP (PF00672), CBS
(PF00571), GGDEF (PF00990), EAL (PF00563), HD-GYP (PF13487). Each entry has
`domain`, `pfam_id`, `chemistry` in `{heme, iron_sulfur, cysteine_rich,
flavin, redox, metal, cnmp, generic_signal}`. Heme and Fe-S domain sets are
hand-curated and committed as part of T0.

### `configs/scoring.yaml`

```yaml
locus_score_weights:
  anchor_marker: 3.0
  accessory_marker: 1.0
  operon_integrity: 1.5
  homology_confidence: 1.0
  taxonomic_context: 0.5
  neighborhood_completeness: 0.5

candidate_score_weights:
  locus: 1.0
  regulator_domain: 1.5
  sensory_domain: 2.0
  proximity: 1.0
  orientation: 0.25
  archetype_conservation: 1.5
  enrichment: 2.0
  taxonomic_breadth: 0.5
  structural_plausibility: 1.5

confidence_thresholds:
  high: 6.0
  medium: 3.5
  low: 1.5

enrichment:
  test: fisher
  multiple_comparison: benjamini_hochberg
  alpha: 0.05
  case_control_ratio: [1, 3]
  permutations: 10000

windows:
  strict: 5
  standard: 10
  extended: 20
  default: standard

robustness:
  windows_to_test: [5, 10, 20]
  weight_perturbation_pct: 20
```

### `configs/benchmarks.yaml`

```yaml
benchmark_csv: data/benchmarks/benchmark_v1.csv
positive_recall_threshold: 0.80
negative_false_positive_threshold: 0.10
report_per_family: true
```

### `configs/headline.yaml`

The single config that drives `make repro`. Composes the analyte configs,
scoring, and benchmark, and pins window=standard, seed=20260429.

## 8. Module specifications

Each subsection lists module purpose, public surface, behavior contract, and
required tests. Function bodies are the agent's responsibility. Function
signatures are not.

### 8.1 `gasregnet/config.py`

Purpose: load and validate YAML configs as Pydantic models.

```python
class AnchorFamilyConfig(BaseModel): ...
class AnalyteConfig(BaseModel): ...
class ScoringConfig(BaseModel): ...
class RegulatorFamilyEntry(BaseModel): ...
class SensoryDomainEntry(BaseModel): ...
class BenchmarkConfig(BaseModel): ...
class GasRegNetConfig(BaseModel):
    analytes: list[AnalyteConfig]
    regulator_families: list[RegulatorFamilyEntry]
    sensory_domains: list[SensoryDomainEntry]
    scoring: ScoringConfig
    benchmark: BenchmarkConfig

def load_config(config_dir: Path) -> GasRegNetConfig: ...
def resolve_and_dump(config: GasRegNetConfig, out_dir: Path) -> Path: ...
```

Tests: round-trip preserves all fields; invalid YAML raises `ConfigError`
naming the offending field path; missing required field raises naming the
field.

### 8.2 `gasregnet/schemas.py`

Purpose: Pandera schemas for every cross-module table.

Exports one schema constant per table from section 6: `LociSchema`,
`GenesSchema`, `RegulatorCandidatesSchema`, `ArchetypesSchema`,
`EnrichmentResultsSchema`, `BenchmarkSchema`. Provides:

```python
def validate(df: pl.DataFrame, schema: pa.DataFrameSchema) -> pl.DataFrame: ...
```

Wraps Pandera and raises `SchemaError` with helpful diagnostics.

Tests: each schema rejects a DataFrame with a missing column, rejects wrong
dtype, accepts a minimal valid DataFrame.

### 8.3 `gasregnet/io/sqlite_efi.py`

Purpose: read EFI-GNT SQLite exports and produce LociSchema, GenesSchema.

```python
def read_efi_sqlite(
    path: Path,
    analyte: str,
    cluster_filter: list[int] | None = None,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Returns (loci_df, genes_df) validated against LociSchema and GenesSchema."""
```

Behavior: reads SSN cluster, neighborhood, and gene tables; emits one `loci`
row per anchor neighborhood and one `genes` row per gene in window. Computes
`relative_index`, `relative_start`, `relative_stop`, `strand`. Sets
`is_anchor=True` for the row matching the SSN query.

Tests: against `tests/fixtures/mini_efi.sqlite`, returns the expected number
of loci and genes; relative indices around the anchor are correct; schema
validation passes.

### 8.4 `gasregnet/io/fasta.py`, `gff.py`, `parquet.py`

Standard helpers. `parquet.py` wraps Polars I/O with schema validation on
read and on write. `gff.py` handles GFF3 only and rejects GFF2. `fasta.py`
yields `(accession, description, sequence)` tuples via BioPython.

Tests: Parquet round-trip preserves schema; GFF reader handles version 3 and
rejects version 2 with a clear error.

### 8.5 `gasregnet/search/diamond.py`

Purpose: run DIAMOND blastp from seeds against a target proteome database.

```python
def run_diamond(
    query_faa: Path,
    db_dmnd: Path,
    out_tsv: Path,
    *,
    evalue: float = 1e-10,
    coverage: float = 0.5,
    identity: float = 0.3,
    threads: int = 8,
    sensitivity: str = "very-sensitive",
) -> Path: ...

def parse_diamond_output(tsv: Path) -> pl.DataFrame: ...
```

Behavior: shells out to the `diamond` binary; captures stderr to logs;
records `diamond --version` in the manifest. Output schema:
`query_id, subject_id, percent_identity, length, mismatch, gapopen, qstart,
qend, sstart, send, evalue, bitscore, qcovhsp, scovhsp`.

Tests: parse output of a known fixture TSV; missing binary raises
`MissingInputError`.

### 8.6 `gasregnet/search/hmmer.py`

```python
def hmmsearch(
    profile_hmm: Path,
    sequences_faa: Path,
    *,
    e_value: float = 1e-5,
) -> pl.DataFrame: ...
```

Implemented through pyhmmer. No subprocess.

### 8.7 `gasregnet/search/mmseqs.py`

```python
def cluster_sequences(
    input_faa: Path,
    out_dir: Path,
    *,
    min_seq_id: float = 0.5,
    coverage: float = 0.8,
) -> pl.DataFrame: ...
```

Returns `sequence_id, cluster_id, is_representative`.

### 8.8 `gasregnet/annotation/domains.py`

```python
def annotate_domains(
    genes: pl.DataFrame,
    pfam_table: pl.DataFrame,
    interpro_table: pl.DataFrame,
) -> pl.DataFrame:
    """Returns genes with pfam_ids, pfam_descriptions, interpro_ids, interpro_descriptions populated."""
```

Idempotent. Logs counts of unannotated genes per locus.

### 8.9 `gasregnet/annotation/regulators.py`

```python
def classify_regulators(
    genes: pl.DataFrame,
    regulator_families: list[RegulatorFamilyEntry],
) -> pl.DataFrame:
    """Sets functional_class, regulator_class, is_regulator_candidate."""
```

Two-component handling is explicit. A gene with `PF00072` (response receiver)
is `two_component_rr`. A gene with `PF00512` (HK domain) or `PF02518`
(HATPase) is `two_component_hk`. Sigma factors (`PF04542`, `PF04545`) are
`sigma`. If no rule fires, `regulator_class = "none"` and
`is_regulator_candidate = False`.

### 8.10 `gasregnet/annotation/taxonomy.py`, `ecology.py`

`taxonomy.py` joins NCBI taxon IDs to phylum, class, order, family, genus
using a local taxonomy dump that the user supplies. `ecology.py` joins to
known-organism tables and produces `taxonomic_context_score` per locus.

### 8.11 `gasregnet/neighborhoods/`

`retrieve.py` produces `loci` and `genes` from a homolog list plus a GFF
directory or an EFI-GNT SQLite. `encode.py` produces an architecture string
per locus following the format
`[-3:HK:PAS][-2:RR:REC-HTH][-1:cydB][0:cydA][+1:cydX]`. `operons.py` infers
operon membership from intergenic distance plus strand co-orientation;
default threshold is 200 nt; configurable.

### 8.12 `gasregnet/scoring/`

`loci.py`: computes `locus_score` from the components defined in
`scoring.yaml`. Each component is a separate function returning a Polars
expression. Final score is the weighted sum.

`candidates.py`: computes `candidate_score` and assigns confidence thresholds
from `scoring.yaml`.

`controls.py`:

```python
def sample_matched_controls(
    case_loci: pl.DataFrame,
    candidate_pool: pl.DataFrame,
    *,
    ratio: tuple[int, int],
    seed: int,
) -> pl.DataFrame: ...
```

Matches on organism and contig, excludes high-confidence anchor loci,
respects `ratio` (default 1:3 case:control).

`enrichment.py`: runs Fisher's exact for each regulator family and sensory
domain feature, then BH-corrects across all features tested per analyte.
Returns `EnrichmentResultsSchema`.

### 8.13 `gasregnet/archetypes/`

`cluster.py`: clusters architecture strings with a distance metric weighted
by proximity to anchor (near genes weighted higher). Use deterministic
agglomerative clustering with a configurable distance threshold. The metric
is documented in the module docstring.

`diagrams.py`: renders gene-arrow figures using matplotlib for each
archetype. Outputs SVG and PNG at 300 DPI. Style file is
`gasregnet/reports/gasregnet.mplstyle`.

### 8.14 `gasregnet/structure/`

`msa.py`: builds MSAs for top candidates using MMseqs2 plus MUSCLE through
BioPython. Returns alignment objects and conserved-residue tables.

`alphafold.py`: ingests AlphaFold-Server outputs (PDB plus pLDDT plus PAE
JSON) and produces a per-residue confidence table. Does not run AlphaFold
itself.

`pdb.py`: aligns AlphaFold models to PDB homologs using a structural
alignment binary if available (TM-align or US-align), with BioPython
superposition as fallback. Produces residue mappings.

### 8.15 `gasregnet/reports/`

`figures.py`: every numbered figure in section 4 of `PREPRINT_PLAN.md` is one
function:

```python
def figure_1_workflow_and_recovery(benchmark_results: pl.DataFrame, out_dir: Path) -> Path: ...
def figure_2_locus_landscape(loci: pl.DataFrame, out_dir: Path) -> Path: ...
def figure_3_archetype_atlas(archetypes: pl.DataFrame, out_dir: Path) -> Path: ...
def figure_4_chemistry_partition(enrichment: pl.DataFrame, out_dir: Path) -> Path: ...
def figure_5_candidate_ranking(candidates: pl.DataFrame, out_dir: Path) -> Path: ...
def figure_6_structure_hypotheses(top_candidates: pl.DataFrame, structures_dir: Path, out_dir: Path) -> Path: ...
```

Outputs SVG and PNG at 300 DPI.

`tables.py`: writes T1 through T6 from section 5 of `PREPRINT_PLAN.md` as
both CSV and Markdown.

`captions.py`: generates a result-led caption sentence per figure, following
the template in section 4 of `PREPRINT_PLAN.md`. Captions are derived from
data, not hard-coded; numbers in captions come from the actual run.

## 9. CLI specification

Typer app. Subcommands:

```bash
gasregnet validate-config --config configs/
gasregnet build-benchmark --out data/benchmarks/benchmark_v1.csv
gasregnet run-sqlite --sqlite path/to/efi.sqlite --analytes CO CN --config configs/ --out results/sqlite_demo
gasregnet diamond-search --query data/seeds/co_anchor_seeds.faa --db databases/bacteria.dmnd --out cache/co_diamond_hits.parquet
gasregnet annotate --neighborhoods <path> --config configs/ --out <path>
gasregnet score --neighborhoods <path> --config configs/ --out <path>
gasregnet enrich --scored <path> --config configs/ --out <path>
gasregnet archetypes --scored <path> --config configs/ --out <path>
gasregnet report --results <dir> --out <dir>
gasregnet repro --config configs/headline.yaml --out results/headline
```

Every subcommand accepts `--verbose`, `--seed`, and writes `manifest.json`
plus `config.resolved.yaml` to its `--out`. `gasregnet --help` lists
everything; each subcommand has a help string.

## 10. Snakemake workflows

`workflows/sqlite_mode.smk` orchestrates the SQLite path. The DAG:

```text
read_sqlite -> annotate_domains -> classify_regulators -> score_loci ->
sample_controls -> enrichment -> score_candidates -> cluster_archetypes ->
figures -> tables
```

Each rule has named inputs and outputs that match the CLI subcommands. Rules
are idempotent and produce Parquet intermediates.

`workflows/diamond_mode.smk` and `workflows/full_discovery.smk` follow the
same downstream graph but replace the head with DIAMOND search plus
neighborhood retrieval from GFFs.

`Makefile` provides:

```text
make sync           # uv sync
make lint           # ruff check + mypy strict
make test           # pytest -q --cov
make repro          # snakemake -s workflows/sqlite_mode.smk on the fixture
make repro-real     # snakemake on the real EFI-GNT SQLite under GNN_analysis_tool/
make clean          # remove results/ and caches
```

## 11. Task graph

Tasks are independently assignable. Definition of done is binding.

### T0. Build the benchmark

Owner: 1 agent. No dependencies.

Output: `data/benchmarks/benchmark_v1.csv` validated against
`BenchmarkSchema`.

Content: at least 30 positive entries split between CO and CN. Required CO
entries include CooA from *Rhodospirillum rubrum* (start from UniProt P19919
class), RcoM-1 and RcoM-2 from *Burkholderia xenovorans*, and at least three
additional literature-curated CO-responsive regulators with PMIDs. Required
CN entries include MpaR from *Pseudomonas aeruginosa* (UniProt A6UZV9),
CydR/Anr-class regulators in *Pseudomonas* with cyanide-tolerant respiration
phenotypes, and CynR-class cyanate regulators distinguished from direct
cyanide sensors in `notes`. Required negative-control entries include
canonical LacI, TetR, AraC, TrpR, LexA from *E. coli* with PMIDs.

DOD: CSV validates against `BenchmarkSchema`; every row has at least one
PMID; CHANGELOG entry; unit test loads the file and confirms minimum counts
per analyte.

### T1. CLAIMS.md

Human-authored. Not a code task. Drives manuscript spine. Three legs of the
thesis, each tied to a figure and a statistical test.

### T2. Repository scaffold

Owner: 1 agent. No dependencies.

Output: directory tree from section 4. `pyproject.toml` with pinned versions
from section 3. `uv.lock`. `Makefile` with stubs. `README.md` with
quickstart pointing to `make repro`. CI workflow at
`.github/workflows/ci.yml`.

DOD: `uv sync` succeeds; `pytest` collects zero tests but does not error;
`ruff check gasregnet/` passes; `mypy --strict gasregnet/` passes on empty
modules.

### T3. Schemas and config models

Owner: 1 agent. Depends on T2.

Output: `gasregnet/schemas.py`, `gasregnet/config.py`, `gasregnet/errors.py`,
`gasregnet/logging.py`, `gasregnet/paths.py`, `gasregnet/hashing.py`,
`gasregnet/manifest.py`. Tests for every schema and every config model.

DOD: every schema in section 6 implemented; every config in section 7 loads,
validates, round-trips; coverage 95% or higher on these modules.

### T4. SQLite ingest

Owner: 1 agent. Depends on T3.

Output: `gasregnet/io/sqlite_efi.py`. `tests/fixtures/mini_efi.sqlite`
produced by `scripts/build_test_fixtures.py`. Integration test ingests the
fixture and validates against `LociSchema`, `GenesSchema`.

DOD: against the EFI-GNT SQLite under `GNN_analysis_tool/`, ingest runs in
under 60 seconds and produces non-empty validated tables.

### T5. Annotation modules

Owner: 1 agent. Depends on T3.

Output: `domains.py`, `regulators.py`, `taxonomy.py`, `ecology.py` complete
with tests.

DOD: applied to T4 output, every gene receives a `functional_class`;
regulator candidates are flagged according to `configs/regulator_families.yaml`;
manual spot-check on top 10 most common Pfam IDs passes.

### T6. Locus scoring

Owner: 1 agent. Depends on T5.

Output: `gasregnet/scoring/loci.py` plus tests.

DOD: deterministic; same input plus same config produces identical scores;
component scores stored as separate columns; weighted total matches the sum
to floating-point tolerance.

### T7. Matched controls and enrichment

Owner: 1 agent. Depends on T6.

Output: `controls.py`, `enrichment.py` plus tests including a deterministic
synthetic case-control fixture where the true odds ratio is known.

DOD: synthetic test recovers the true odds ratio within 5% across 100 random
seeds; BH correction matches scipy and statsmodels reference values.

### T8. Candidate scoring

Owner: 1 agent. Depends on T7.

Output: `candidates.py` plus tests. Confidence-class assignment using
thresholds from `scoring.yaml`.

DOD: ranking output reproduces a hand-computed top-5 from a synthetic locus
set; every score-decomposition column is populated; rationale strings are
non-empty.

### T9. Archetype clustering and diagrams

Owner: 1 agent. Depends on T8.

Output: `archetypes/cluster.py`, `archetypes/diagrams.py` plus tests.

DOD: clustering deterministic given the seed; gene-arrow diagrams render as
SVG and PNG for each archetype with the architecture string embedded as a
caption.

### T10. Benchmark recovery integration test

Owner: 1 agent. Depends on T0, T4, T5, T6, T8.

Output: `tests/integration/test_benchmark_recovery.py`. Runs the SQLite
pipeline on `GNN_analysis_tool/`'s SQLite and reports recall on the benchmark
positive set and false-positive rate on the negative set. Writes
`results/benchmark/recovery.csv`.

DOD: integration test exists, runs in under 5 minutes, asserts recall above
the threshold in `configs/benchmarks.yaml`.

### T11. Partition claim integration test

Owner: 1 agent. Depends on T7, T8.

Output: `tests/integration/test_partition_claim.py`. Runs enrichment on CO
and CN separately, computes whether sensory chemistries differ between
analytes via chi-square on the joint family-by-chemistry table at q < 0.05.

DOD: test runs to completion; outcome is logged whether the partition holds;
result written to `results/headline/partition_outcome.json`. Failure of the
partition does not fail the test; it surfaces a finding.

### T12. Reports and figures

Owner: 1 agent. Depends on T8, T9, T10, T11.

Output: every figure function in section 8.15 produces a 300 DPI PNG and an
SVG with a result-led caption. Every table in section 5 of
`PREPRINT_PLAN.md` produces both CSV and Markdown.

DOD: `make repro` produces `results/headline/figures/` containing six PNGs
and six SVGs and `results/headline/tables/` containing T1 through T6.

### T13. Snakemake workflow

Owner: 1 agent. Depends on T2 through T9, T12.

Output: `workflows/sqlite_mode.smk` complete. `make repro` invokes Snakemake
with the SQLite fixture and produces the headline figures and tables.

DOD: clean clone plus `uv sync` plus `make repro` reproduces the headline
outputs end to end in under 30 minutes on a developer laptop.

### T14. Robustness panel

Owner: 1 agent. Depends on T8, T7.

Output: a script that runs the pipeline at window sizes 5, 10, 20 and at
+/- 20% weight perturbations; outputs Spearman concordance of the top-100
candidate list across perturbations.

DOD: concordance table written to `results/headline/robustness.csv`;
concordance below 0.7 emits a warning in the manifest.

### T15. Structure module skeleton

Owner: 1 agent. Depends on T8.

Output: `structure/msa.py`, `alphafold.py`, `pdb.py` with implementations
covering ingestion of AlphaFold-Server outputs from a fixture and structural
alignment to a fixture PDB.

DOD: ingestion of fixture AlphaFold output produces a per-residue confidence
table; alignment to fixture PDB produces a residue mapping; no online calls
made.

### T16. CLI

Owner: 1 agent. Depends on T2, T3, T4, T6, T8, T9, T13.

Output: `gasregnet/cli.py` complete. Every subcommand from section 9
implemented and tested with `typer.testing.CliRunner`.

DOD: `gasregnet --help` lists all subcommands; each subcommand has a help
string and at least one CLI test.

### T17. Manifest and hashing

Owner: 1 agent. Depends on T2, T3.

Output: `gasregnet/hashing.py`, `gasregnet/manifest.py`. Every CLI subcommand
calls them.

DOD: every results directory contains `manifest.json` with documented schema;
running with the same inputs twice produces the same input and config
hashes.

### T18. Documentation

Owner: 1 agent. Depends on T13, T16.

Output: `README.md` with quickstart; `docs/methods.md` describing the methods
text appropriate for the manuscript; `docs/reproducibility.md` describing
the Zenodo deposit plan.

DOD: a researcher unfamiliar with the repo can run `make repro` and produce
headline outputs from the README alone.

## 12. Testing requirements

Three tiers.

Unit tests live in `tests/unit/<module>/`. Coverage on each module marked P1
or P2 in the task graph must be 80% or higher. Use Hypothesis for schema
validation, scoring monotonicity, and BH correction monotonicity.

Integration tests live in `tests/integration/`. Required: SQLite end-to-end
on a fixture; benchmark recovery; partition claim. Each integration test
produces a results directory with `manifest.json`.

End-to-end is `make repro` itself. CI must run `make repro` on the fixture
and assert on output file presence and on a checksum of the headline figures
table.

CI is GitHub Actions, two jobs. Job 1: `uv sync`, `ruff check`, `mypy --strict`,
`pytest -q --cov`. Job 2: `make repro` on the fixture. Job 1 is required for
merge. Job 2 is informational at first and required after T13.

## 13. Anti-patterns

Do not introduce a new DataFrame library. Polars is the library. Pandas is
allowed only at integration boundaries.

Do not write to disk outside of `--out`. Test fixtures are read-only.

Do not silently catch exceptions. Catch typed, log with context, re-raise or
convert to a typed exception and raise.

Do not make schemas optional. Validate at every read and every write.

Do not introduce a deep-learning dependency without explicit go-ahead.
Heuristic, statistical, and structural approaches are sufficient.

Do not chase coverage with tautological tests. A test that asserts `True` is
worse than no test.

Do not bypass Snakemake by chaining CLI calls in `Makefile` rules other than
`make repro`.

Do not rely on internet access except in modes that explicitly opt in.

Do not edit `ARCHITECTURE.md`, `PREPRINT_PLAN.md`, or `ENGINEERING_PROMPT.md`
as part of a code task. If a contract is wrong, raise it as an issue and let
the human revise.

Do not extend the inherited scripts in `GNN_analysis_tool/`. They are
historical context. New code lives under `gasregnet/`.

## 14. Definition of done at the repository level

The repository is preprint-ready when all of the following hold.

1. `make repro` runs end to end in under 30 minutes on a developer laptop
   using the SQLite fixture and produces the six headline figures and six
   headline tables.
2. The benchmark recovery integration test reports recall above the
   configured threshold and false-positive rate below it.
3. The partition claim integration test has been run and the outcome
   recorded in `results/headline/partition_outcome.json`. If the partition
   does not hold, the human revises the manuscript thesis. The codebase
   remains DOD-met.
4. Every Pandera schema in section 6 is enforced at every cross-module
   boundary.
5. Coverage on `gasregnet/` is 80% or higher; mypy strict is clean; ruff is
   clean.
6. Every figure has a result-led caption committed alongside it.
7. The repository has a tagged release at `v0.1.0` with a Zenodo DOI and a
   `CITATION.cff` file.
8. `manifest.json` accompanying the headline run records the package version,
   all config hashes, all input data hashes, the seed, the wall clock, and
   the version of every external binary used.

When all eight conditions are met, the codebase is the preprint's Methods
section made executable.

## 15. Parallelization plan

Five-agent slice after T0 through T3 land:

- Agent A: T2, then T13 (scaffold, then Snakemake)
- Agent B: T0, then T10 (benchmark, then benchmark recovery test)
- Agent C: T4, T5 (ingest and annotation)
- Agent D: T6, T7, T8 (scoring chain)
- Agent E: T9, T12 (archetypes and reports)

After this slice:

- Agent A: T16 (CLI), T17 (manifest)
- Agent B: T11 (partition test), T14 (robustness)
- Agent C: T15 (structure)
- Agent D: T18 (docs)
- Agent E: support and integration glue

T3 is the gating task. No agent makes meaningful progress without the
schemas. Land T3 before launching the wider fan-out.
