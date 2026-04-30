# GasRegNet Architecture

GasRegNet is a full-scale comparative genomics platform for discovering candidate
gas-responsive bacterial transcriptional regulators. The inherited GNN scripts are
treated as evidence of the original scientific goal, not as the architectural
boundary of the new project.

The central hypothesis is:

> Bacterial genes required for CO metabolism or cyanide tolerance are often
> embedded in conserved local genomic architectures. Regulatory proteins that
> repeatedly co-occur near those loci, carry plausible sensory domains, and are
> enriched relative to matched controls are strong candidate gas sensors.

GasRegNet should therefore infer:

```text
gas-associated functional locus
  -> conserved genome neighborhood
  -> regulator candidate
  -> ranked sensor hypothesis
  -> structure-guided experimental target
```

## Current Inherited Assets

The current repository contains:

- `GNN_analysis_tool/gnn_sqlite_analyzer.py`
  - Reads an EFI-GNT SQLite export.
  - Scores CO dehydrogenase-like neighborhoods using cox marker Pfams.
  - Detects nearby transcription-factor-like Pfams.
  - Produces cluster summaries, TF tables, archetype patterns, and figures.
- `GNN_analysis_tool/gnn_analyzer.py`
  - Older TSV-oriented analyzer.
  - Useful mostly as historical context.
- `GNN_analysis_tool/config/pfam_config.json`
  - CO marker and broad TF-family substring rules.
- `GNN_analysis_tool/config/known_co_oxidizers.csv`
  - Literature-backed CO oxidizer/tolerator reference table.
- `CHM6635_FinalProject_Info_Rubric_ML_addendum.pdf`
  - Defines the class deliverable, but not the ceiling of this project.

These assets establish Dr. Dent's scientific direction: use functional genes as
anchors, inspect neighboring genes, identify putative regulators, and rank
sensor candidates.

## First-Principles Scope

GasRegNet is not a CO-only GNN post-processor. It is an analyte-general discovery
system. The first supported analytes are:

- `CO`: carbon monoxide metabolism / aerobic CO oxidation.
- `CN`: cyanide tolerance, cyanide-resistant respiration, and cyanide detox.

The architecture should also support future analytes such as NO, H2S, O2,
methane, ammonia, or other small-molecule gas signals by adding configuration and
seed data rather than rewriting the pipeline.

## Scientific Questions

GasRegNet should answer:

1. Which loci are high-confidence CO- or CN-associated functional loci?
2. Which SSN or homology clusters contain those loci?
3. Which local gene architectures recur across genomes?
4. Which regulator families are enriched near high-confidence gas loci?
5. Which individual regulatory proteins are the strongest sensor candidates?
6. Do CO and CN loci show distinct regulatory architectures?
7. Which candidates have structural features consistent with direct gas sensing?
8. What specific experiments should test the top-ranked sensor hypotheses?

## System Overview

```text
seed proteins / HMM profiles / EFI-GNT SQLite / proteome databases
        |
        v
homology search and anchor detection
        |
        v
genome neighborhood retrieval
        |
        v
domain, regulator, taxonomy, and ecology annotation
        |
        v
locus scoring and confidence classification
        |
        v
operon archetype encoding and clustering
        |
        v
matched-control enrichment statistics
        |
        v
candidate regulator ranking
        |
        v
publication-style figures, tables, and report assets
        |
        v
structure-guided shortlist and validation proposals
```

## Evidence Layers

Each candidate sensor is scored from independent evidence layers. Scores must be
decomposable so the final ranking is interpretable.

### 1. Functional Anchor Layer

Identify gas-associated functional genes.

CO anchor examples:

- `coxL`: aerobic carbon monoxide dehydrogenase large subunit.
- `coxM`: medium/FAD subunit.
- `coxS`: small Fe-S subunit.
- `coxG`, `coxI`, and other accessory/chaperone markers when supported.

CN anchor examples:

- cyanide-tolerant terminal oxidase markers.
- cytochrome bd oxidase components such as `cydA`, `cydB`, `cydX`.
- cyanide detox or tolerance-associated genes when biologically justified.
- sulfurtransferase/rhodanese, nitrilase/cyanide hydratase, cyanase/cyanate
  pathway genes as contextual markers, not automatic proof.

The anchor layer should distinguish:

- direct functional marker evidence.
- accessory/contextual marker evidence.
- weak homology evidence.
- ambiguous oxidase evidence.

### 2. Homology Layer

Use scalable protein search to expand beyond manually exported EFI-GNT data.

Primary tools:

- DIAMOND for large-scale protein similarity search.
- HMMER for marker/profile validation.
- MMseqs2 for clustering and redundancy reduction.

The homology layer should support:

- seed FASTA searches.
- reciprocal-best-hit checks.
- profile-HMM confirmation.
- cluster assignment.
- coverage, identity, e-value, and bitscore filtering.

### 3. Genome Context Layer

For each anchor hit, retrieve local neighborhoods.

Required fields:

- anchor accession.
- contig/genome identifier.
- organism and taxonomic metadata.
- neighboring gene accessions.
- relative gene position.
- nucleotide distance from anchor.
- strand/orientation.
- predicted operon membership where possible.
- neighborhood boundary/truncation flags.

Neighborhood windows should be configurable:

```text
strict window: +/- 5 genes
standard window: +/- 10 genes
extended window: +/- 20 genes
```

### 4. Regulatory Candidate Layer

Detect nearby regulatory genes.

Regulator categories:

- one-component transcription factors.
- two-component response regulators.
- sensory histidine kinases.
- sigma and anti-sigma systems.
- transcriptional antiterminators.
- riboswitches or RNA regulators as a future module.

Important regulator families:

- FNR/CRP.
- LysR.
- GntR.
- AraC.
- MarR.
- TetR.
- MerR.
- LacI.
- ArsR/SmtB.
- LuxR/GerE.
- XRE.
- response regulators.
- histidine kinases.

### 5. Sensory Domain Layer

Classify domains that make a regulator biologically plausible as a gas sensor.

Examples:

- PAS.
- GAF.
- Cache.
- HAMP.
- CBS.
- heme-binding domains.
- flavin/redox domains.
- Fe-S-associated domains.
- metal-binding domains.
- GGDEF/EAL/HD-GYP cyclic-di-GMP domains.
- receiver domains in two-component systems.

For CO and CN, special attention should be paid to domains or motifs involving:

- heme.
- iron-sulfur centers.
- cysteine/histidine metal coordination.
- flavin/redox chemistry.
- oxygen/redox sensing.
- electrophile/nucleophile-sensitive residues.

### 6. Ecology And Taxonomy Layer

Add organism-level context.

CO evidence:

- known carboxydotrophs.
- known carboxydovores.
- literature-backed CO oxidizers.

CN evidence:

- known cyanide producers.
- known cyanide-tolerant organisms.
- organisms with cyanide-resistant respiration.
- organisms from cyanide-rich ecological contexts where supported.

Taxonomy must be used carefully to avoid over-weighting clonal or genus-specific
signals. Enrichment tests should include taxonomically matched controls whenever
possible.

### 7. Conservation Layer

Measure recurrence across genomes and homolog clusters.

Important features:

- frequency of regulator co-occurrence near high-confidence loci.
- conservation of regulator class within an archetype.
- conservation across taxonomic breadth.
- proximity to anchor gene.
- strand/orientation conservation.
- specificity to CO or CN loci relative to controls.

### 8. Structure Layer

For top candidates only:

- build or retrieve AlphaFold/ColabFold models.
- map DNA-binding and sensory domains.
- construct MSA for homologous candidate regulators.
- map conserved residues onto the structure.
- identify PDB homologs.
- align structures.
- nominate residues or pockets plausibly involved in gas sensing.

Structure is not a bulk step. It is a final prioritization layer for the top
ranked candidates.

## Scoring Model

Each locus receives a gas-associated functional score. Each regulator receives a
candidate sensor score.

### Locus Score

```text
locus_score =
  anchor_marker_score
+ accessory_marker_score
+ operon_integrity_score
+ homology_confidence_score
+ taxonomic_context_score
+ neighborhood_completeness_score
```

Locus confidence classes:

- `high`: multiple marker classes, coherent local architecture, strong homology.
- `medium`: convincing anchor with partial supporting context.
- `low`: weak or ambiguous evidence.
- `control`: non-gas or low-confidence comparison locus.

### Regulator Candidate Score

```text
candidate_score =
  locus_score
+ regulator_domain_score
+ sensory_domain_score
+ proximity_score
+ orientation_score
+ archetype_conservation_score
+ enrichment_score
+ taxonomic_breadth_score
+ structural_plausibility_score
```

Scores must be stored as separate columns as well as a final weighted total.
This prevents the final ranking from becoming a black box.

## Operon Archetype Discovery

Each neighborhood should be encoded as a centered gene architecture:

```text
[-5] [TF:LysR] [-4] [unknown] [-3] [coxM] [-2] [coxS] [-1] [coxG] [0] [coxL] [+1] [coxI]
```

or:

```text
[-3] [HK:PAS] [-2] [RR:REC-HTH] [-1] [cydB] [0] [cydA] [+1] [cydX] [+2] [sulfurtransferase]
```

Similarity should be weighted by distance from the anchor:

- near genes matter most.
- regulator class and sensory domain matches matter strongly.
- exact accession identity matters less than domain/function class.
- unknown genes contribute weakly.
- strand/orientation similarity adds support.

Outputs:

- `archetypes.csv`
- `archetype_members.csv`
- representative locus for each archetype.
- frequency by analyte and cluster.
- taxonomic breadth.
- gene-arrow diagrams for top archetypes.

## Matched Controls

Controls are required for strong claims.

Control types:

1. Matched non-gas oxidase loci.
2. Low-confidence anchor hits.
3. Random genes from the same genomes.
4. Regulators far from the anchor, such as 20-50 genes away.
5. Taxonomically matched background loci.

Questions:

- Are specific regulator families enriched near high-confidence CO loci?
- Are specific regulator families enriched near high-confidence CN loci?
- Are sensory domains enriched near high-confidence gas loci?
- Are regulator distances shorter than expected by chance?
- Are archetypes more conserved than random neighborhoods?

Statistical methods:

- Fisher's exact test.
- odds ratios.
- permutation tests.
- false-discovery-rate correction.
- taxonomic stratification where feasible.

## Data Model

### `loci`

One row per anchor neighborhood.

Required columns:

```text
locus_id
analyte
anchor_accession
anchor_family
organism
taxon_id
cluster_id
contig_id
window_size
is_boundary_truncated
marker_genes_present
accessory_genes_present
locus_score
locus_confidence
taxonomic_context_score
operon_integrity_score
```

### `genes`

One row per gene in a neighborhood.

Required columns:

```text
locus_id
gene_accession
relative_index
relative_start
relative_stop
strand
product_description
pfam_ids
pfam_descriptions
interpro_ids
interpro_descriptions
functional_class
regulator_class
sensory_domains
is_anchor
is_regulator_candidate
```

### `regulator_candidates`

One row per putative regulatory protein near a gas-associated locus.

Required columns:

```text
candidate_id
analyte
locus_id
gene_accession
organism
cluster_id
relative_index
distance_nt
position
strand
regulator_class
dna_binding_domains
sensory_domains
pfam_ids
interpro_ids
archetype_id
locus_score
regulator_domain_score
sensory_domain_score
proximity_score
archetype_conservation_score
enrichment_score
taxonomic_breadth_score
structural_plausibility_score
candidate_score
rationale
```

### `archetypes`

One row per recurring architecture.

Required columns:

```text
archetype_id
analyte
cluster_id
architecture_string
n_loci
n_taxa
representative_locus_id
dominant_regulator_class
dominant_anchor_structure
mean_locus_score
mean_candidate_score
```

### `enrichment_results`

One row per statistical comparison.

Required columns:

```text
analyte
feature_type
feature_name
case_definition
control_definition
n_case_with_feature
n_case_without_feature
n_control_with_feature
n_control_without_feature
odds_ratio
p_value
q_value
interpretation
```

## Proposed Repository Structure

```text
GasRegNet/
  ARCHITECTURE.md
  README.md
  pyproject.toml
  configs/
    analytes/
      co.yaml
      cn.yaml
    regulator_families.yaml
    sensory_domains.yaml
    scoring.yaml
  data/
    seeds/
      co_anchor_seeds.faa
      cn_anchor_seeds.faa
    references/
      known_co_organisms.csv
      known_cn_organisms.csv
    controls/
  gasregnet/
    cli.py
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
  workflows/
    sqlite_mode.smk
    diamond_mode.smk
    full_discovery.smk
  results/
    loci.csv
    genes.csv
    regulator_candidates.csv
    archetypes.csv
    enrichment_results.csv
    figures/
    report_assets/
```

## CLI Design

SQLite demonstration mode:

```bash
gasregnet run-sqlite \
  --sqlite GNN_analysis_tool/*.sqlite \
  --analytes CO CN \
  --out results/sqlite_demo
```

Full discovery mode:

```bash
gasregnet discover \
  --analytes CO CN \
  --proteome-db databases/bacteria.dmnd \
  --genome-metadata metadata/genomes.tsv \
  --gff-dir genomes/gff \
  --out results/full_discovery
```

DIAMOND-only expansion:

```bash
gasregnet diamond-search \
  --query data/seeds/co_anchor_seeds.faa \
  --db databases/bacteria.dmnd \
  --out cache/co_diamond_hits.parquet
```

Figure/report generation:

```bash
gasregnet report \
  --results results/full_discovery \
  --out results/full_discovery/report_assets
```

## Execution Modes

### Mode 1: EFI-GNT SQLite Mode

Purpose:

- reproduce and greatly extend Dr. Dent's GNN-style analysis.
- produce immediate class-project deliverables.
- validate scoring, regulator classification, and figure generation.

Inputs:

- EFI-GNT SQLite file.
- analyte configs.
- marker/domain configs.
- known organism tables.

Outputs:

- scored loci.
- regulator candidates.
- archetypes.
- enrichment tables where controls can be defined from the SQLite data.
- publication-style figures.

### Mode 2: Curated Proteome Panel Mode

Purpose:

- scalable but bounded analysis.
- useful for debugging the full pipeline.
- includes known positives, known negatives, and taxonomic controls.

Inputs:

- curated proteome FASTA.
- genome annotation files.
- seed proteins.

Outputs:

- DIAMOND hits.
- retrieved neighborhoods.
- same scored outputs as SQLite mode.

### Mode 3: Full DIAMOND Discovery Mode

Purpose:

- real-scale homolog discovery.
- preprint-oriented expansion beyond the EFI-GNT export.

Inputs:

- large DIAMOND database.
- genome metadata.
- genome annotation files or a resolvable protein-to-genome coordinate map.
- seed proteins and HMM profiles.

Outputs:

- full homolog universe.
- clustered anchor families.
- scored neighborhoods.
- regulator enrichment and candidate rankings.

## Major Figures

### Figure 1: GasRegNet Workflow

Pipeline schematic from seed proteins and genome data to ranked gas-sensor
candidates.

### Figure 2: CO And CN Locus Landscape

Cluster-level landscape showing high-, medium-, and low-confidence loci for CO
and CN. The query cluster should be marked when an SSN is used.

### Figure 3: Archetypal Gene Architectures

Gene-arrow diagrams for top CO and CN archetypes, with frequency and taxonomic
breadth.

### Figure 4: Regulator Enrichment

Heatmap or dot plot showing regulator families and sensory domains enriched near
high-confidence CO and CN loci.

### Figure 5: Candidate Sensor Ranking

Top candidates with score decomposition: locus, proximity, regulator domain,
sensory domain, enrichment, conservation, and structure.

### Figure 6: Structure-Guided Sensor Hypotheses

One top CO candidate and one top CN candidate. Show AlphaFold/PDB comparison,
domain architecture, conserved residues, and proposed sensing pocket.

## Tables

Required high-impact tables:

- major clusters and high-confidence loci.
- all putative regulators across major clusters.
- regulator classification by Pfam, InterPro, class, and sensory domain.
- top ranked candidate sensors.
- enrichment statistics.
- representative archetypes.

## Report / Manuscript Frame

The project should be written as:

> We developed GasRegNet, a modular comparative genomics pipeline for discovery
> of gas-responsive bacterial regulatory proteins. GasRegNet integrates homology
> search, genome-neighborhood architecture, domain annotation, operon archetype
> clustering, matched-control enrichment, and structure-guided prioritization.
> Applied to CO and cyanide-associated loci, the platform identifies distinct
> regulatory architectures and nominates specific transcription factors as
> experimentally testable candidate gas sensors.

## Implementation Priorities

Priority 1:

- define analyte configs for CO and CN.
- produce normalized `loci`, `genes`, and `regulator_candidates` tables.
- score loci and candidates.
- classify regulators and sensory domains.

Priority 2:

- implement archetype encoding and clustering.
- generate gene-arrow diagrams.
- produce all rubric-relevant tables and figures.

Priority 3:

- add controls and enrichment statistics.
- produce CO-vs-CN comparative figures.
- write automated captions and methods text.

Priority 4:

- add DIAMOND-backed full discovery mode.
- add HMMER validation.
- add structure-guided shortlist assets.

## Design Principles

- Config-driven analyte support.
- Reproducible outputs from raw inputs.
- Interpretable scores, not opaque black-box predictions.
- Controls before claims.
- Store intermediate data in durable tabular formats.
- Keep SQLite mode as a fast test fixture.
- Make full DIAMOND discovery optional but real.
- Prioritize ranked biological hypotheses over raw hit lists.

