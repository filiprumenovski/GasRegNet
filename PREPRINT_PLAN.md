# GasRegNet to bioRxiv: Preprint Reframing Plan

A complement to `ARCHITECTURE.md`. This document does not replace the
architecture. It changes the contract the architecture is being built to
satisfy. The architecture currently builds a tool. A preprint requires a
finding. Every section below is opinionated and intended to be acted on.

## 1. The transition is not a polish job

ARCHITECTURE.md answers "how is this built." A preprint must answer "what was
discovered, why is it true, why does anyone care." Those are different
documents. The current spec frames GasRegNet as a discovery framework that
produces ranked candidates. That framing yields a methods paper, but a weak
one, because reviewers see no falsifiable result, only a tool. A *Bioinformatics*
or *PLOS Computational Biology* methods reviewer will accept it. A reviewer at
any biological venue will not. The pipeline becomes Methods. The findings
become Results. The discovery becomes the abstract. Every figure caption,
table title, and supplement description has to be rewritten so it terminates
in a biological sentence rather than a capability statement.

## 2. The thesis you should commit to

> Bacterial sensors of CO and HCN partition into chemically stratified
> regulatory architectures. Systematic genome-neighborhood mining recovers
> known sensors, expands each family by N-fold, and nominates a structurally
> tractable shortlist of previously uncharacterized candidates.

Three independently falsifiable legs:

1. **Recovery.** Recall of canonical bacterial gas-responsive regulators is
   above a stated threshold on a versioned benchmark set. CooA, RcoM, FNR/CRP
   redox-sensing variants, NorR-family, DosS/DosT, MpaR, CydR/Anr-class
   regulators belong on this benchmark. Negative controls (LacI, TetR
   canonical, AraC canonical) belong on it too.

2. **Discovery.** Across a defined corpus, the method nominates novel
   candidates whose ranking is enriched for sensory domains and conserved
   metal-coordinating residues relative to matched non-gas oxidase loci at
   FDR < 0.05.

3. **Insight.** CO-anchored loci enrich for heme and Fe-S coordination within
   FNR/CRP and CooA-like clades. Cyanide-anchored loci enrich for cysteine-rich
   MerR/ArsR families and two-component systems flanking cytochrome bd. The
   partition is the headline biological result.

Without leg 3, the paper is a tool report and will be cited as such. With leg
3, the paper carries a finding that travels.

## 3. Title and abstract

Title: **Genome-neighborhood mining partitions bacterial CO and HCN sensors by
coordination chemistry.**

Abstract, target ~225 words:

> Bacteria sense small toxic gases including carbon monoxide and hydrogen
> cyanide using transcriptional regulators that have been characterized one
> family at a time. The full landscape of gas-responsive regulators across
> bacterial diversity is unknown, and there is no systematic resource that
> ranks plausible sensor candidates with calibrated confidence. We developed
> GasRegNet, a comparative genomics framework that integrates homology search,
> genome-neighborhood encoding, sensory-domain priors, operon archetype
> clustering, and matched-control enrichment to nominate bacterial gas-sensing
> regulators. Applied to N anchor families spanning CO oxidation and
> cyanide-tolerant respiration across M sequenced bacterial genomes from GTDB
> rs220, GasRegNet recovers all canonical CO sensors at recall X and identifies
> P high-confidence regulator candidates passing matched-control enrichment at
> q < 0.05. CO-anchored and cyanide-anchored loci partition into chemically
> distinct regulatory architectures. CO loci concentrate heme and iron-sulfur
> sensory chemistry within FNR/CRP and CooA-like clades, while cyanide-anchored
> loci enrich for cysteine-rich MerR/ArsR families and two-component systems
> flanking cytochrome bd. AlphaFold-guided structural prioritization nominates
> one CO and one cyanide candidate with conserved metal-coordinating residues
> lining a closed pocket positioned to bind a diatomic ligand. We release
> ranked candidate tables, archetype catalogs, structural models, and a
> reproducible Snakemake pipeline as a community resource.

Treat the abstract as a contract. Every Methods choice you have not yet made
is forced by it.

## 4. The five gaps ARCHITECTURE.md does not cover

### Gap 1. Benchmark and positive control set

There is no named benchmark in the current spec. Without a list of known true
sensors that the pipeline must recover, no reviewer can assess whether the
ranking is meaningful. Build a versioned benchmark of 30 to 60 literature-curated
bacterial gas-responsive regulators with PMID, organism, UniProt accession,
anchor gene, and direct or indirect sensing claim. Include matched negative
regulators of comparable family abundance that are well-characterized but not
gas-responsive. Report recall on positives, false-positive rate on negatives,
and per-family sensitivity. Without this asset, no candidate ranking is
publishable. This is the highest-leverage thing to produce next, before
anything else.

### Gap 2. Statistical backbone

Current spec mentions Fisher's exact and FDR correction without committing.
Reviewers will require, and you should commit to:

- A named null model. State explicitly that for each regulator family the null
  is occurrence within ±10 genes of high-confidence anchor loci at the same
  rate as occurrence within ±10 genes of taxonomically matched non-gas oxidase
  loci.
- A named control corpus, sized 1:3 case to control, sampled to match
  organism-level abundance.
- Benjamini-Hochberg correction across all regulator families at q < 0.05,
  with the family count stated.
- A robustness panel showing rankings are stable across ±5, ±10, and ±20 gene
  windows, reported as Spearman concordance of the top-100 candidate list.
- Permutation tests for archetype enrichment at 10,000 permutations minimum.

ARCHITECTURE.md schedules statistics in Priority 3. For a preprint, statistics
is load-bearing from Priority 1. Move it.

### Gap 3. Comparative positioning

GasRegNet is not the first tool to mine genome neighborhoods. Reviewers will
demand explicit positioning against EFI-GNT and EFI-EST as the proximate
prior, RODEO for operon-context discovery, antiSMASH as the reference for
preprint-scale cluster mining, FlaGs and webFlaGs for flanking-gene
visualization, gggenes and gggenomes for the visualization stack, and
DeepBGC-class tools for ML-driven cluster scoring.

Build a feature-comparison table for the manuscript with tools on the x-axis
and features on the y-axis: analyte-general configuration, decomposable
scoring, matched-control enrichment, archetype clustering, structural
prioritization, FDR-controlled candidates, reproducible workflow. The honest
differentiation is the integration of decomposable scoring with matched-control
enrichment and structural prioritization in a single configurable workflow.
That is real. State it explicitly. Bury it and reviewers will dismiss the
paper as "EFI-GNT plus a wrapper."

### Gap 4. Figures must terminate in findings, not capabilities

ARCHITECTURE.md lists six figures, all of which describe what the pipeline
shows rather than what was discovered. Reframe each caption to lead with a
result sentence:

- **Figure 1.** *GasRegNet recovers known bacterial CO and HCN sensors across
  N organisms with recall X.* Workflow schematic in panel A, ROC and PR
  curves on the benchmark in panels B and C, score-component contribution
  bars in panel D.

- **Figure 2.** *CO and HCN anchor loci occupy distinct taxonomic territory
  and cluster into separable architectures.* Locus landscape colored by
  analyte and confidence class with the query cluster marked.

- **Figure 3.** *Six recurrent archetypes account for Y% of high-confidence
  CO loci and six for Z% of HCN loci.* Gene-arrow atlas with frequencies
  and taxonomic breadth.

- **Figure 4.** *Sensory chemistry partitions CO and HCN regulator families.*
  Regulator family by sensory domain heatmap, faceted by analyte, with
  enrichment q values overlaid. This is the headline figure.

- **Figure 5.** *GasRegNet nominates P high-confidence candidate sensors with
  decomposable scores.* Top-30 ranking with stacked-bar score decomposition.

- **Figure 6.** *Structural prioritization identifies a closed metal-coordinated
  pocket in the top CO candidate and a cysteine-lined pocket in the top HCN
  candidate.* AlphaFold model, PDB superposition, conserved-residue map,
  proposed binding pocket.

If a caption only describes the visualization, the figure has not earned its
slot.

### Gap 5. Limitations and falsifiability

Pre-empt the four standard objections in a half-page Limitations subsection:

1. **Annotation propagation error.** Many regulator-family assignments come
   from automated Pfam and InterPro propagation. Report agreement across
   Pfam, InterPro, and HMMER profile hits. Flag candidates that depend on a
   single annotation source.

2. **Co-occurrence is not regulation.** Physical proximity does not prove
   transcriptional control. State this directly. GasRegNet generates
   regulator-target *hypotheses*, not *evidence*. The class rubric's
   Biochemistry experiment becomes Discussion in the preprint, recast as a
   prospective validation program.

3. **Genome corpus bias.** GTDB and RefSeq overweight clinically and
   industrially relevant taxa. Report phylum-level coverage. Discuss what
   may be missing from undersampled phyla.

4. **CO and HCN are not the same biology.** Lumping them as "small gas
   molecules" risks obscuring real differences. Lead with the partition,
   not the unification.

A clean Limitations section before Discussion materially improves review
outcomes.

## 5. Tables, reorganized for preprint use

| Table | Content |
|-------|---------|
| T1 | Benchmark recovery: known-sensor list, hit/miss, rank, score decomposition |
| T2 | Top 30 CO candidates with anchor, organism, regulator class, sensory domain, score components, q value, structural-pocket flag |
| T3 | Top 30 HCN candidates, same schema |
| T4 | Regulator-family enrichment by analyte with FDR |
| T5 | Archetype catalog: frequency, taxonomic breadth, dominant regulator class |
| T6 | Tool feature-comparison matrix |

All per-locus and per-gene tables move to supplements. T1, T2, T3, T6 are
non-negotiable in the main text.

## 6. Methods rigor required for a preprint

ARCHITECTURE.md mentions Snakemake, Parquet, configs. For preprint-grade
reproducibility, also commit to:

- A Zenodo deposit at submission containing the pinned conda or mamba
  lockfile, the Snakemake DAG export, configs used for the headline run,
  input data hashes, and the exact GTDB or RefSeq snapshot tag.
- A `make repro` target that reproduces all main figures from a single
  command on a clean clone in under a defined wall-clock budget.
- Versioned releases starting at v0.1.0 for the preprint version, with
  CITATION.cff and a Zenodo DOI cited in the preprint.
- A Methods subsection naming every external tool with version: DIAMOND,
  HMMER, MMseqs2, InterProScan, AlphaFold-Server snapshot date, PyMOL,
  Cytoscape if used.

This is non-negotiable for a methods-heavy preprint. Reviewers will check.

## 7. Implementation priority resequencing

ARCHITECTURE.md priorities are correct for a class deliverable. For a preprint,
resequence around what most reduces preprint risk:

| Phase | Window | Deliverable |
|-------|--------|-------------|
| P0 | week 0 | Benchmark positive and negative gene set, versioned CSV in `data/benchmarks/` |
| P1 | week 1 to 2 | Locus and gene tables, deterministic decomposable scoring, positive-control recall as a CI test |
| P2 | week 2 to 3 | Matched-control sampling, Fisher and BH, robustness panels, ROC and PR on benchmark |
| P3 | week 3 to 4 | Archetype clustering, gene-arrow figures, chemistry-partition figure with q values |
| P4 | week 4 to 5 | Top-candidate ranking, AlphaFold for top 2 to 4 per analyte, PDB alignment, PyMOL figures |
| P5 | week 5 to 6 | Manuscript draft, methods rigor pass, supplements, clean-clone reproducibility, Zenodo deposit |

DIAMOND full-discovery mode stays below the preprint waterline unless the
EFI-GNT corpus lacks statistical power for the chemistry-partition claim. Run
a power calculation early in P0. If GTDB scale is required, push DIAMOND into
P2.

## 8. Risk register

Four risks most likely to kill the preprint in review, ranked by likelihood:

1. **The chemistry partition does not survive matched controls.** Run the
   partition test in P2 before figure design. If it fails, the paper becomes
   a tools paper at the *Bioinformatics* level, which is publishable but
   requires reframing. Better to learn this in week 3 than week 6.

2. **Known-positive recall is poor because the EFI-GNT SQLite corpus does
   not contain enough canonical sensor neighborhoods.** Validate benchmark
   coverage at the end of P0. If CooA neighborhoods are not in the SQLite
   export, switch to DIAMOND discovery before scoring.

3. **AlphaFold structural claims for top candidates are non-distinctive.** Low
   pTM or high PAE in the binding region kills Figure 6. Pre-compute pLDDT
   and PAE for the top 10 in each analyte during P3, not P4. Choose
   figure-eligible candidates only after structural QC.

4. **Reviewer dismissal as "EFI-GNT plus a wrapper."** Front-load
   differentiation in Introduction and the comparison table. Make matched-control
   enrichment, decomposable scoring, and analyte-general configuration the
   explicit novelty statements rather than buried features.

## 9. What to write next, in order

1. The benchmark positive and negative gene set as a versioned CSV in
   `data/benchmarks/`. Highest-leverage asset to produce immediately.
2. A one-page `CLAIMS.md` listing the three legs of the thesis, each with
   the figure that supports it and the statistical test that justifies it.
   This becomes the spine of the manuscript.
3. A `power_calc.ipynb` that estimates whether the SQLite corpus has the
   statistical power to detect the chemistry partition at the target FDR. If
   yes, run the preprint in SQLite mode. If no, schedule DIAMOND.
4. A first-pass abstract committed to the repo so every subsequent design
   choice can be tested against "does this support the abstract."

The architecture is sound. The discovery has not yet been made. Make the
discovery first. The manuscript follows from it cleanly.
