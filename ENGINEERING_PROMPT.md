# GasRegNet V2 Engineering Prompt: Scientific Audit and Correction Plan

A layered addition to `ENGINEERING_PROMPT.md`. V1 stands. V2 corrects the
biology, raises the analytical layer to publication standard, and
introduces the modules required to clear bioRxiv-grade peer review. Read
`ARCHITECTURE.md`, `PREPRINT_PLAN.md`, and `ENGINEERING_PROMPT.md` first.
Then read this in full before writing code.

This document is opinionated. Every prescription has a biological reason.
Where the reason is non-obvious, the rationale is given inline. Agents
implement what is specified, not what they think would be reasonable.

## 0. How to use this document

Sections 1 and 2 are the audit. Read them so you understand why V2 exists.
Section 3 is required biology reading; the seed citations there are how the
benchmark and configs were constructed and how the agents must reason about
edge cases. Section 4 is the V2 mission delta. Sections 5 and 6 are pinned
configurations. Section 7 is module-by-module specification. Sections 8
through 13 are contracts, tasks, and definition-of-done.

V1 is not deprecated. V1 schemas remain valid. V2 extends them with new
optional columns and adds new tables. Where V2 corrects a V1 default, the
correction is explicit in section 6.

If anything in this document conflicts with V1, V2 wins for the file or
contract under discussion. If V2 leaves something unspecified, V1 wins. No
agent reinvents a contract on its own.

## 1. Audit summary: what V1 got right

Before the corrections, credit where due. The engineering layer is strong
and should not be rewritten. V2 builds on top.

The V1 codebase enforces typed Pandera schemas at every cross-module
boundary. Polars is used end-to-end with no pandas leakage. DuckDB is
already used as the genome catalog backend. Snakemake DAGs orchestrate the
SQLite and corpus modes. Manifest hashing, configuration validation through
Pydantic, structlog logging, mypy strict mode, and per-module tests are all
present. The repository can be cloned, synced through `uv`, and exercised
end-to-end through `make repro`. This is real infrastructure and V2 will
keep all of it.

The schema for `loci`, `genes`, `regulator_candidates`, `archetypes`,
`enrichment_results`, and `benchmark` is correct in shape. The locus
scoring decomposes cleanly. The partition test is wired to a real chi-square
on an analyte by domain table. Architecture-string clustering is
deterministic.

What V1 lacks is biological correctness. That is the entire scope of V2.

## 2. Audit summary: fourteen scientific defects in V1

These are the defects V2 must fix. Each is named, located, and given the
required correction. Numbered for cross-reference from later sections.

### Defect 1. The benchmark is anchor genes, not regulators

`data/benchmarks/benchmark_v1.csv` lists `cooS`, `coxL`, `coxM`, `coxS`,
`cydA`, `cydB`, `cydX` as benchmark entries. These are the *target genes*
that gas-responsive regulators control. They are not the regulators
themselves. Recovery of an anchor gene that the term-scan already searches
for is tautological. The pipeline cannot fail this benchmark because the
benchmark is a self-test of the scanner.

Required correction: rebuild the benchmark with literature-curated bacterial
gas-responsive regulators. The canonical list is in section 6.1. Every
entry has a UniProt accession, a published characterization, and, where
applicable, a structural or biochemical demonstration of direct gas
sensing.

### Defect 2. Anchor detection is a substring match on product strings

`gasregnet.datasets.refseq.scan_refseq_catalog` runs
`product ILIKE '%coxL%'` against the catalog. The branch that requires
`mode='smoke'` raises `NotImplementedError` for anything else. This is
text matching on annotation strings, not anchor detection. Pseudomonas
genomes annotated with `cytochrome bd-II ubiquinol oxidase` are matched
by the `cytochrome bd` pattern but represent different cyanide-tolerance
biology than canonical cydAB cytochrome bd-I. Genomes with synonyms or
abbreviations are missed.

Required correction: profile-driven anchor detection. HMMER profiles
hit against translated CDS sequences with bitscore and per-family
gathering thresholds set by the curator. DIAMOND back-search against
seed proteins for homology confirmation. Term-scan retained only as a
permissive recall-augmenting layer with a separate
`evidence_type='term_match'` flag. See section 7.1.

### Defect 3. Sensory chemistry assignments are biologically wrong

`configs/sensory_domains.yaml` maps every PAS variant to `heme`. PAS
(PF00989, PF13426, PF14598) is a structural fold that holds many ligand
chemistries. Among PAS domains in bacteria, FAD-binding (NifL, Aer),
4-hydroxycinnamate-binding (PYP family), and citrate-binding cases are at
least as common as heme-binding. Heme-PAS is the minority case and is
specifically the FixL, EcDOS, RcoM-N family. Tagging all PAS as heme
inflates `sensory_domain_score` for every PAS-bearing protein.

GAF (PF01590) is mapped to `generic_signal`. GAF is canonically a
cyclic-nucleotide binding domain. The DosS, DosT, and NreB GAF domains
bind heme; the cGMP phosphodiesterase GAFs bind cGMP. The correct
chemistries are `cnmp` and `heme`, not `generic_signal`.

HAMP (PF00672) is a coiled-coil signal transduction module that lies
between sensor and output domains. It does not bind ligand. Listing
HAMP as a sensory domain inflates scores for transduction proteins
that do not sense.

CBS (PF00571) is mapped to `redox`. CBS pairs bind adenosine derivatives
(AMP, ATP, S-adenosylmethionine). They are allosteric regulators, not
redox sensors.

GGDEF (PF00990), EAL (PF00563), and HD-GYP (PF13487) are output domains
for c-di-GMP synthesis and hydrolysis. They are response effectors, not
sensors. They are misclassified as `generic_signal` in V1, which both
inflates the sensor score and confuses sensors with effectors.

H-NOX (PF07700, also called HNOB), the canonical heme-NO/oxygen binding
domain that defines the bacterial nitric oxide and CO sensor family in
many genomes, is absent from the config.

Globin (PF00042) and the globin-coupled sensor architecture (HemAT and
relatives) are absent from the config. These sense O2 and bind CO and NO
through the canonical heme-iron mechanism.

Required correction: replace `configs/sensory_domains.yaml` with the
chemistry-stratified, ligand-validated mapping in section 6.2. Add a
distinct enumeration `effector_domain` for GGDEF, EAL, HD-GYP, and
remove them from sensor scoring.

### Defect 4. Sensors and regulators are conflated under one role

V1 marks any gene matching a Pfam in `regulator_families.yaml` as a
`regulator_candidate` and any gene with a sensory-domain Pfam contributes
`sensory_domain_score` to that candidate. The biological reality has three
roles, not one. A gene can be:

- An anchor: the gas-utilizing functional gene at position 0.
- A regulator: a DNA-binding transcription factor in the locus.
- A sensor: the protein that physically binds the gas ligand. The sensor
  may be the regulator itself (CooA, RcoM, FNR), the cognate histidine
  kinase of a two-component regulator (FixL of FixJ, DosS of DosR), or a
  separate accessory sensor in the operon.

V1 treats two-component HKs as regulators (`regulator_class =
two_component_hk`) when in fact the HK is the sensor and the response
regulator is the DNA binder. This makes the partition claim untestable
because the chemistry of the sensor is collapsed onto the partition of
the regulator.

Required correction: introduce a `sensor_role` column on `genes`
(`anchor`, `regulator`, `sensor`, `accessory`, `none`) independent of
`regulator_class`. Score sensor evidence on the sensor row. Score
regulator evidence on the regulator row. Pair them through a new
`sensor_regulator_pairs` table when both are present in the same locus.
See sections 6.4 and 7.2.

### Defect 5. Regulator classification is first-match, not evidence-accumulating

`gasregnet.annotation.regulators._classify_regulator` returns on first hit:
PF00072 -> `two_component_rr`, then PF00512 or PF02518 -> `two_component_hk`,
then PF04542 or PF04545 -> `sigma`, then loops over families. A FNR/CRP
homolog with a fragmentary HK annotation from automated propagation gets
classified as `two_component_hk` and never reaches the FNR/CRP rule. A
sigma-54 enhancer-binding protein with a receiver domain is classified as
`two_component_rr` and its NorR-class identity is lost.

Required correction: evidence-accumulating classifier. Compute a
multi-class score vector from all matched Pfams. Resolve ties using
priority: family-specific HTH-bearing rules win over generic receiver
or kinase rules unless no other DNA-binding domain is present. See
section 7.2.

### Defect 6. PDB structural alignment is index-zip, not alignment

`gasregnet.structure.pdb.residue_mapping_by_order` zips two CA atom lists
by sequence order. For homologs of equal length and identical alignment,
this works by coincidence. For homologs of any non-trivial divergence, the
returned mapping is meaningless. Any conserved-residue identification on
top of this mapping is also meaningless. Figure 6 is currently
unpublishable on this implementation.

Required correction: real structural alignment via TM-align or US-align.
Both are fast, free, and produce a residue-residue mapping with TM-score
and aligned RMSD. Wrap the binary in `gasregnet.structure.align`. Where
the binary is unavailable, fall back to a sequence-then-structure
pipeline using a Biopython superposer applied to residues paired by an
MSA-derived alignment. Never zip CA atoms by order. See section 7.6.

### Defect 7. operon_integrity_score and taxonomic_context_score are zero everywhere

In `gasregnet.io.sqlite_efi` and `gasregnet.datasets.refseq`, every
emitted locus row sets `operon_integrity_score = 0.0` and
`taxonomic_context_score = 0.0`. The locus scorer multiplies by weights
1.5 and 0.5 respectively and the result contributes nothing. The real
operon integrity computation in `gasregnet.neighborhoods.operons` is not
wired into the locus output. The taxonomic context computation in
`gasregnet.annotation.ecology` is not invoked.

Required correction: wire the existing `operons.anchor_operon_integrity`
and `ecology.score_taxonomic_context` outputs into the locus emit step in
both ingest paths. Where ecology data is absent for an organism, emit
`null` and let the scorer treat null as zero with a manifest warning.

### Defect 8. Conservation layer is unimplemented

V1 sets `archetype_conservation_score = 0.0` for every candidate.
Conservation is the strongest evidence for a real regulatory architecture
and is named as Layer 7 in `ARCHITECTURE.md`. Without it, archetype
clustering exists but contributes nothing to the candidate ranking, and
the headline preprint claim that recurrent architectures predict gas
sensing is unsupportable.

Required correction: compute `archetype_conservation_score` as a function
of the count of distinct genera the candidate's archetype spans, the
fraction of high-confidence loci in the analyte cluster bearing the same
architecture token at the same relative index, and a minimum-loci
threshold below which the score is zero. See section 7.4.

### Defect 9. Enrichment has no phylogenetic correction

`gasregnet.scoring.enrichment.run_enrichment` runs Fisher's exact test on
case versus control loci treating each locus as one independent
observation. Twenty Pseudomonas loci from one clade contribute as twenty
independent observations of FNR/CRP enrichment. Reviewers will reject
this. The correct statistic stratifies by genus or family and combines
either through Cochran-Mantel-Haenszel or through phyla-level
deduplication.

Required correction: stratified enrichment by configurable taxonomic
rank (default genus). Cochran-Mantel-Haenszel for the per-feature
combined odds ratio with continuity correction. Robustness panel runs
the same enrichment under three deduplication policies (none, one
representative per genus, one representative per family) and reports
concordance. See section 7.5.

### Defect 10. Distance is not in nucleotides

`gasregnet.scoring.candidates._distance_nt` uses
`min(abs(relative_start), abs(relative_stop))` where `relative_start` and
`relative_stop` come from V1's neighborhood ingest. In `sqlite_efi` these
are gene-index distances, not nucleotide coordinates. The column name and
the proximity score that consumes it are misleading. Any biological
inference that uses absolute distances or operon co-transcription
thresholds is wrong.

Required correction: maintain two columns. `relative_index` is the gene
index (already correct). `distance_nt` is the absolute nucleotide
distance from the anchor gene's stop to the candidate gene's start, or
from the candidate's stop to the anchor's start, whichever is smaller
and on the same strand. When neither is available because contig
coordinates are missing, set `distance_nt = null` and propagate the
null through proximity scoring with a manifest warning.

### Defect 11. Architecture-string distance over-weights anchor identity

`gasregnet.archetypes.cluster.architecture_distance` weights position 0
at full weight (1.0). The anchor token at position 0 then dominates any
distance between two architectures that differ only in anchor family.
Two cydA loci with identical flanking architecture group together; a
cydA and a cydB locus with otherwise identical architecture do not. This
collapses cytochrome bd-I and bd-II loci into separate archetypes when
their regulatory biology is closely related.

Required correction: separate anchor-family distance from flanking-gene
distance. When two architectures share a curator-defined anchor family
group, anchor-position cost is zero. Curator-defined groups are in
section 6.3 (for example `cyd_oxidase_group = {cydA, cydB, appC, appB}`).
For across-group comparisons, anchor mismatch contributes weight 1.0 as
in V1.

### Defect 12. cydX Pfam is wrong and cydX is structurally hard to detect

V1 sets `cydX` to `pfam_required: ['PF13456']`. PF13456 is not the
Pfam family for CydX. CydX is a 37-residue accessory subunit whose
Pfam coverage is incomplete in many bacterial genomes. Term-scan
on `cydX` will miss it in any genome where the locus tag and product
fields do not contain the string. The accessory marker presence count
will be artificially low.

Required correction: add a fallback rule. A short open reading frame
(less than 60 residues) co-localized within 200 nucleotides of a
cydA-cydB pair on the same strand is annotated as `cydX_like` and
contributes to accessory marker count with a lower confidence weight.
Pfam ID is changed to PF14392 if the curator confirms the assignment;
otherwise the field is left empty and the inference is rule-based. See
section 6.3.

### Defect 13. Seed FASTAs mix aerobic and anaerobic CO dehydrogenase families

`data/seeds/co_anchor_seeds.faa` contains P19919 (Afipia coxL, aerobic
Mo-CODH large subunit) and P31896 (Rhodospirillum rubrum cooS, anaerobic
Ni-Fe-CODH catalytic subunit). These are different enzyme families with
different cofactors, different Pfam membership, and different evolutionary
ancestry. A DIAMOND search seeded with both retrieves a non-monophyletic
set. The downstream cluster assignment treats them as one family.

Required correction: split CO seeds into `co_aerobic_seeds.faa` and
`co_anaerobic_seeds.faa`. Run anchor detection separately for each.
Maintain analyte `CO` as the union but tag each anchor hit with
`anchor_subfamily in {aerobic_codh, anaerobic_codh}`. See section 6.3.

### Defect 14. Operator-region motif evidence is absent

A regulator candidate adjacent to an anchor is a co-occurrence
hypothesis. A regulator candidate with a recognizable binding site in
the anchor's operator region is a regulation hypothesis. V1 generates
only the former. Without operator evidence, the preprint cannot
distinguish "this regulator is near this gene" from "this regulator
likely controls this gene". This is the central biological claim of the
paper and it has no module backing it.

Required correction: add a motif and operator-region module. Extract
the upstream intergenic region of each anchor (default minus 200 to
plus 50 relative to the anchor start). Run de novo motif discovery
through STREME on the high-confidence locus set per analyte. Match
discovered motifs against curated TF binding sites in CollecTF and
RegPrecise. Score per-locus operator support and add it to candidate
scoring. See section 7.3.

## 3. Required reading and authoritative biology references

These are the papers V2 was constructed from. The benchmark and config
files cite them. Before any agent edits a config, the agent reads at least
the abstract of the relevant references.

CO sensing in bacteria. Aono et al. 1996 J. Biol. Chem. 271:31963-31968,
discovery of CooA in *Rhodospirillum rubrum*. Lanzilotta et al. 2000 Nat.
Struct. Biol. 7:876-880, CooA crystal structure with heme. Roberts et al.
2004 Microbiol. Mol. Biol. Rev. 68:453-473, review of CO sensing through
heme proteins. Marvin et al. 2008 J. Biol. Chem. 283:23194-23201,
characterization of RcoM-1 and RcoM-2 from *Burkholderia xenovorans*
LB400. King and Weber 2007 Nat. Rev. Microbiol. 5:107-118, ecology of
aerobic CO oxidation.

Aerobic CO dehydrogenase. Dobbek et al. 2002 Proc. Natl. Acad. Sci. U.S.A.
99:15971-15976, structure of *Oligotropha carboxidovorans* CO
dehydrogenase. The aerobic Mo-Cu CODH is a CoxLMS heterotrimer; coxL
contains the active site. Anaerobic Ni-Fe CODH is the cooS family and is
unrelated in fold and cofactor.

Cyanide tolerance in bacteria. Cunningham et al. 1997 Microbiology
143:3873-3880, cytochrome bd as the cyanide-resistant terminal oxidase in
*Pseudomonas aeruginosa*. Borisov et al. 2011 Biochim. Biophys. Acta
1807:1398-1413, review of cytochrome bd biochemistry. Pessi and Haas
2000 J. Bacteriol. 182:6940-6949 and Pessi and Haas 2001 FEMS
Microbiol. Lett. 200:73-78, regulation of HCN biosynthesis in
*P. aeruginosa* through ANR. Castric 1994 Curr. Microbiol. 29:19-21
and Blumer and Haas 2000 Arch. Microbiol. 173:170-177, cyanogenesis
genetics. Carterson et al. 2004 J. Bacteriol. 186:6837-6844, MpaR-class
regulators.

Two-component and one-component sensors. Galperin 2018 Annu. Rev.
Microbiol. 72:347-368, taxonomy of bacterial signal transduction. Ulrich
et al. 2005 Trends Microbiol. 13:52-56, one-component versus
two-component systems. Mascher et al. 2006 Microbiol. Mol. Biol. Rev.
70:910-938, sensor histidine kinases.

Heme-PAS and heme-GAF sensors. Gilles-Gonzalez and Gonzalez 2005 J.
Inorg. Biochem. 99:1-22, heme-based sensors. Sasakura et al. 2002 J.
Biol. Chem. 277:23821-23827, EcDOS as a heme-PAS oxygen sensor. Sousa
et al. 2007 Protein Sci. 16:1708-1719, DosS/DosT GAF heme binding in
*M. tuberculosis*. Kumar et al. 2007 Proc. Natl. Acad. Sci. U.S.A.
104:11568-11573, DosT senses oxygen.

H-NOX domain. Iyer et al. 2003 BMC Genomics 4:5, identification of
HNOB / H-NOX as a bacterial heme-NO/O2 sensor family. Pellicena et al.
2004 Proc. Natl. Acad. Sci. U.S.A. 101:12854-12859, *Thermoanaerobacter*
H-NOX structure.

Pfam, InterPro, and HMMER. El-Gebali et al. 2019 Nucleic Acids Res.
47:D427-D432, Pfam 32.0 (and the rolling Pfam release thereafter).
Eddy 2011 PLoS Comput. Biol. 7:e1002195, HMMER methodology.

DIAMOND. Buchfink et al. 2021 Nat. Methods 18:366-368, DIAMOND v2 with
sensitive and ultra-sensitive modes.

ESM-2. Lin et al. 2023 Science 379:1123-1130, evolutionary-scale protein
language modeling and the ESM-2 family.

Foldseek. van Kempen et al. 2024 Nat. Biotechnol. 42:243-246,
ultra-fast structure-based protein search.

TM-align and US-align. Zhang and Skolnick 2005 Nucleic Acids Res.
33:2302-2309 (TM-align). Zhang et al. 2022 Nat. Methods 19:1109-1115
(US-align).

Operator and motif resources. Bailey et al. 2015 Nucleic Acids Res.
43:W39-W49, MEME suite. Bailey 2021 Bioinformatics 37:2834-2840,
STREME for de novo motif discovery. Kilic et al. 2014 Nucleic Acids
Res. 42:D156-D160, CollecTF database. Novichkov et al. 2013 BMC
Genomics 14:745, RegPrecise 3.0.

Statistics. Mantel and Haenszel 1959 J. Natl. Cancer Inst. 22:719-748,
the test that bears their name. Benjamini and Hochberg 1995 J. R. Stat.
Soc. B 57:289-300, FDR. Felsenstein 1985 Am. Nat. 125:1-15, on the
phylogenetic non-independence problem.

When in doubt, the agent searches Pfam and InterPro for the current
accession of a domain rather than relying on the IDs hard-coded here.
Pfam version drift over the past five years has been minor but
non-zero, and any one-off correction should be flagged in the V2
manifest.

## 4. V2 mission delta

V1 mission: build a tool that ranks gas-sensor candidates with
decomposable scores. V2 mission: build a tool whose ranked candidates
are biologically credible enough to commit one or two for experimental
validation in a wet-lab program. The bar moves from "the pipeline
runs and produces a list" to "a domain expert reading the list would
spend a postdoc-year on the top three".

Concretely V2 requires:

1. A literature-curated regulator benchmark with at least thirty
   positive entries, ten per analyte plus negative controls.
2. Profile-based anchor detection through HMMER and DIAMOND.
3. Three-role gene annotation (anchor, regulator, sensor).
4. Sensory-chemistry classification corrected to match canonical
   ligand binding.
5. Operon integrity, taxonomic context, and archetype conservation
   scores actually computed and contributing to candidate ranking.
6. Phylogenetically corrected enrichment.
7. Operator-region motif discovery and binding-site matching.
8. Real structural alignment with a real algorithm.
9. ESM-2 embedding similarity as an independent evidence channel.
10. Foldseek structural homology against the AlphaFold database.
11. Bayesian scoring with weights fit on the benchmark by
    cross-validation.
12. DuckDB analytical kernel exposing every candidate, score
    component, archetype, and enrichment row as a queryable view.

Items 1 through 8 are P0 and gate the preprint. Items 9 through 12
are P1 and elevate the preprint from a tools paper to a method
paper with a finding.

## 5. Tech stack delta

Add to V1's pinned stack. All additions go through `uv` and
`pyproject.toml`. External binaries are pinned by version recorded in
the manifest at runtime, not at install time.

```text
DataFrame:                       polars 1.x        (V1; unchanged)
Schema validation:               pandera 0.23.x    (V1; unchanged)
Workflow:                        snakemake 8.x     (V1; unchanged)
Storage:                         duckdb 1.x        (V1; unchanged but role expanded)
HMM:                             pyhmmer 0.10.x    (V1; unchanged but newly invoked)
Statistics:                      scipy, statsmodels (V1; unchanged)

NEW Python:
Embeddings:                      fair-esm 2.0.0    (ESM-2; CPU and GPU)
                                 torch 2.4.x      (CPU build by default)
Motif discovery:                 pymeme 0.x is unstable; use the meme binary
                                 with a thin subprocess wrapper. Parse outputs.
Foldseek bindings:               foldseek binary 8.ef4e960 (external)
TM-align:                        tmalign / USalign binary 20240319 (external)
Bayesian fitting:                pymc 5.x for the headline model
                                 scikit-learn 1.5.x for cross-validated
                                 logistic regression as a fast baseline
Provenance:                      prov 2.x for PROV-O metadata export

NEW external binaries:
HMMER:                           hmmer 3.4         (system or conda)
DIAMOND:                         diamond 2.1.10    (V1 listed; now required)
MMseqs2:                         mmseqs 15-6f452   (V1 listed; now required)
MEME suite:                      meme 5.5.7
Foldseek:                        foldseek 8.ef4e960
US-align:                        USalign 20240319
InterProScan:                    interproscan 5.x  (optional for genome-scale)
```

The optional `gpu` extra in `pyproject.toml` brings in
`torch-cuda`-equivalent wheels for ESM-2 acceleration. Default
installation runs CPU-only and ESM-2 large is gated behind the GPU
extra.

Add a `tools.yaml` config file that records the path and version of
every external binary. The CLI exposes a `gasregnet check-tools`
subcommand that probes each binary, records `--version`, and writes
`results/<run>/tools_resolved.yaml` at the start of every run.

DuckDB role expansion. V2 promotes DuckDB from "storage backend for
RefSeq catalogs" to "analytical kernel". The reasons are concrete: the
joins in enrichment and archetype clustering are cleaner in SQL than in
Polars iter_rows; biologists who cannot read Polars can read SQL; and
DuckDB queries Parquet files in place without ingest. The Polars layer
remains for in-memory transformation. DuckDB views are the queryable
interface for biologists, downstream notebooks, and the web UI in V2.

## 6. Corrected configurations

Each subsection is a full file replacement. Where a V1 file is replaced,
the V1 file is kept under `configs/legacy/v1/` for diff comparison.

### 6.1 `configs/benchmarks/regulators_v2.csv`

Replaces `data/benchmarks/benchmark_v1.csv` for any test that purports to
measure regulator recovery. The V1 file is kept for anchor recovery.

Schema follows `BenchmarkSchema` from V1 (`gasregnet.schemas`). At
minimum thirty positive entries are required: ten CO regulators, ten CN
regulators, ten negative controls. The agent commits the file with an
explicit `verify_pmid` column flagged where the PMID was not directly
confirmed and must be checked before manuscript submission.

Required CO regulators (analyte=CO):

| protein_name | uniprot_accession | organism | anchor_family | expected_regulator_class | expected_sensory_domains | sensing_evidence_class | first_publication |
|---|---|---|---|---|---|---|---|
| CooA | P72324 | Rhodospirillum rubrum | cooS | one_component | heme | direct | Aono 1996 JBC |
| RcoM-1 | Q13EP1 | Burkholderia xenovorans LB400 | coxL | one_component | heme_pas | direct | Marvin 2008 JBC |
| RcoM-2 | Q13EP2 | Burkholderia xenovorans LB400 | coxL | one_component | heme_pas | direct | Marvin 2008 JBC |
| FixL | P23222 | Bradyrhizobium diazoefficiens | none_in_locus | two_component_hk | heme_pas | direct | Gilles-Gonzalez 1991 Nature |
| EcDOS | P76129 | Escherichia coli K-12 | none_in_locus | one_component | heme_pas | direct | Sasakura 2002 JBC |
| DosS | P9WGK1 | Mycobacterium tuberculosis | dosR_regulon | two_component_hk | heme_gaf | direct | Sousa 2007 Protein Sci |
| DosT | P9WGK3 | Mycobacterium tuberculosis | dosR_regulon | two_component_hk | heme_gaf | direct | Kumar 2007 PNAS |
| HemAT | O30739 | Bacillus subtilis | none_in_locus | one_component | heme_globin | direct | Hou 2000 PNAS |
| NreB | Q7DLR9 | Staphylococcus carnosus | nitrate_reductase | two_component_hk | iron_sulfur_pas | direct | Mullner 2008 J Mol Microbiol Biotechnol |
| FNR | P0A9E5 | Escherichia coli K-12 | various | one_component | iron_sulfur_4Fe4S | direct | Khoroshilova 1997 PNAS |

Required CN regulators (analyte=CN):

| protein_name | uniprot_accession | organism | anchor_family | expected_regulator_class | expected_sensory_domains | sensing_evidence_class | first_publication |
|---|---|---|---|---|---|---|---|
| ANR | P25856 | Pseudomonas aeruginosa PAO1 | various incl cyd | one_component | iron_sulfur_4Fe4S | direct | Sawers 1991 Mol Microbiol |
| MpaR-class regulator | A6UZV9 | Pseudomonas aeruginosa PA7 | mhp operon | one_component | cysteine_metal | indirect | Microbiology 2010 (rubric ref) |
| ArcA | P0A9Q1 | Escherichia coli K-12 | cydAB | two_component_rr | none | indirect | Iuchi 1988 PNAS |
| ArcB | P0AEC3 | Escherichia coli K-12 | cydAB | two_component_hk | quinone | indirect | Iuchi 1990 J Bacteriol |
| Fnr (cyanide context) | P0A9E5 | Escherichia coli K-12 | cydAB regulation | one_component | iron_sulfur_4Fe4S | indirect | Cotter 1990 J Bacteriol |
| RoxR | Q88FX6 | Pseudomonas putida KT2440 | terminal oxidases | two_component_rr | none | indirect | Comolli 2002 J Bacteriol |
| RoxS | Q88FX5 | Pseudomonas putida KT2440 | terminal oxidases | two_component_hk | redox_pas | indirect | Comolli 2002 J Bacteriol |
| CydR (= ANR ortholog) | (Azotobacter) | Azotobacter vinelandii | cydAB | one_component | iron_sulfur_4Fe4S | direct | Wu 2000 Mol Microbiol |
| CynR | P27111 | Escherichia coli K-12 | cynTSX | one_component | none | indirect | Anderson 1990 J Bacteriol |
| GcsR | (Pseudomonas) | Pseudomonas aeruginosa | hcn operon | sigma54 | aaa_atpase | indirect | Carterson 2004 J Bacteriol |

Required negative controls (analyte=negative_control):

| protein_name | uniprot_accession | organism | expected_regulator_class | rationale |
|---|---|---|---|---|
| LacI | P03023 | Escherichia coli K-12 | one_component | classical lactose regulator; no gas biology |
| TetR | P04483 | Escherichia coli (transposon) | one_component | tetracycline regulator |
| AraC | P0A9E0 | Escherichia coli K-12 | one_component | arabinose regulator |
| TrpR | P0A881 | Escherichia coli K-12 | one_component | tryptophan regulator |
| LexA | P0A7C2 | Escherichia coli K-12 | one_component | SOS response, not gas |
| GalR | P03024 | Escherichia coli K-12 | one_component | galactose regulator |
| RpoH | P0AGB3 | Escherichia coli K-12 | sigma | heat shock sigma factor |
| BglG | P11546 | Escherichia coli K-12 | antiterminator | sugar antitermination |
| GlpR | P0ACL0 | Escherichia coli K-12 | one_component | glycerol regulator |
| MalT | P06993 | Escherichia coli K-12 | one_component | maltose regulator |

Each row has `pmid` filled with the canonical first publication where
known, and `verify_pmid` set when the agent did not directly confirm. The
agent commits both columns.

The benchmark CSV is written by `gasregnet build-benchmark --version v2`
which composes this table from `configs/benchmarks/regulators_v2.csv` and
validates against `BenchmarkSchema`. The V1 anchor benchmark stays at
`data/benchmarks/anchors_v1.csv` and is renamed in place.

### 6.2 `configs/sensory_domains.yaml` (full replacement)

```yaml
# V2 sensory and cofactor domain catalog.
# chemistry values: heme, iron_sulfur_4Fe4S, iron_sulfur_2Fe2S,
#                   flavin, globin, cnmp, cysteine_metal, redox_quinone,
#                   metal_zinc, none
# role values: sensor, transducer, effector
# A domain entry combines a Pfam ID with the curator's chemistry call.
# Chemistry assignments are based on the dominant ligand observed in
# characterized representatives, not on the bare structural fold.

sensory_domains:
  - domain: PAS_generic
    pfam_id: PF00989
    role: sensor
    chemistry: none
    notes: "Generic PAS fold. Default chemistry is none. Heme-PAS, FAD-PAS, and ligand-PAS are scored only when paired with co-occurring chemistry evidence (heme-binding histidine motif, FAD/FMN cofactor signature, etc.)."
  - domain: PAS_9
    pfam_id: PF13426
    role: sensor
    chemistry: none
    notes: "PAS subfamily; default none unless paired evidence."
  - domain: PAS_11
    pfam_id: PF14598
    role: sensor
    chemistry: none
    notes: "PAS subfamily."
  - domain: GAF
    pfam_id: PF01590
    role: sensor
    chemistry: cnmp
    notes: "Default chemistry cyclic-nucleotide; heme-GAF (DosS/DosT/NreB lineage) is rescored to heme by paired evidence rule."
  - domain: HNOB
    pfam_id: PF07700
    role: sensor
    chemistry: heme
    notes: "H-NOX heme NO/O2 binding domain. Iyer 2003 BMC Genomics."
  - domain: Globin
    pfam_id: PF00042
    role: sensor
    chemistry: heme
    notes: "Globin fold; HemAT-like globin-coupled sensors."
  - domain: Cache_1
    pfam_id: PF02743
    role: sensor
    chemistry: none
    notes: "Periplasmic small-molecule sensor; default chemistry none."
  - domain: Cache_2
    pfam_id: PF13185
    role: sensor
    chemistry: none
  - domain: HAMP
    pfam_id: PF00672
    role: transducer
    chemistry: none
    notes: "Coiled-coil signal transduction; not a sensor. Excluded from sensor scoring."
  - domain: CBS
    pfam_id: PF00571
    role: sensor
    chemistry: none
    notes: "Adenosine-derivative binding; not a redox sensor."
  - domain: GGDEF
    pfam_id: PF00990
    role: effector
    chemistry: none
    notes: "Diguanylate cyclase output domain; excluded from sensor scoring."
  - domain: EAL
    pfam_id: PF00563
    role: effector
    chemistry: none
  - domain: HD-GYP
    pfam_id: PF13487
    role: effector
    chemistry: none
  - domain: Fer4
    pfam_id: PF00037
    role: sensor
    chemistry: iron_sulfur_4Fe4S
    notes: "Generic 4Fe-4S binding; FNR-class iron-sulfur sensor."
  - domain: Fer4_7
    pfam_id: PF12838
    role: sensor
    chemistry: iron_sulfur_4Fe4S
  - domain: Fer2
    pfam_id: PF00111
    role: sensor
    chemistry: iron_sulfur_2Fe2S
  - domain: Fer2_2
    pfam_id: PF01799
    role: sensor
    chemistry: iron_sulfur_2Fe2S
    notes: "Common 2Fe-2S; downweighted because it is highly generic."
  - domain: FAD_binding
    pfam_id: PF00890
    role: sensor
    chemistry: flavin
  - domain: cNMP_binding
    pfam_id: PF00027
    role: sensor
    chemistry: cnmp
    notes: "CooA contains this domain but binds heme, not cAMP. Heme rescore applied through CooA-class rule (see paired evidence below)."
  - domain: ArsR_HTH
    pfam_id: PF01022
    role: sensor
    chemistry: cysteine_metal
    notes: "ArsR/SmtB metalloregulator family; senses metal ions through cysteine clusters."

paired_evidence:
  - rule_name: cooA_heme_rescore
    if_pfam_all: ["PF00027"]
    if_motif_any: ["heme_b_histidine"]
    rescore:
      domain: cNMP_binding
      chemistry: heme
    notes: "CooA-class. cNMP fold with proximal His ligation rescores to heme."
  - rule_name: dos_gaf_heme
    if_pfam_all: ["PF01590"]
    if_co_pfam_any: ["PF00512", "PF02518"]
    if_organism_kingdom: bacteria
    rescore:
      domain: GAF
      chemistry: heme
    notes: "GAF in HK context biased toward heme-GAF (DosS, DosT, NreB lineage)."
```

The pipeline applies sensor-domain assignments first, then runs the
`paired_evidence` rules to rescore in light of co-occurring evidence.
Effector and transducer domains never contribute to `sensory_domain_score`.

### 6.3 `configs/analytes/co.yaml` and `configs/analytes/cn.yaml` (corrected)

`co.yaml`:

```yaml
analyte: CO
display_name: "carbon monoxide"
anchor_seeds:
  aerobic_codh: data/seeds/co_aerobic_seeds.faa
  anaerobic_codh: data/seeds/co_anaerobic_seeds.faa
anchor_subfamilies:
  - name: aerobic_codh
    seeds: aerobic_codh
    primary_pfam: PF02738   # Mo-CODH large subunit dimerization
    supporting_pfam:
      - PF03450             # CO dehydrogenase flavoprotein middle subunit
      - PF01799             # Fer2_2 small subunit (downweighted, generic)
    hmm_profiles:
      - data/profiles/coxL.hmm
      - data/profiles/coxM.hmm
      - data/profiles/coxS.hmm
    seed_uniprot: [P19919]
  - name: anaerobic_codh
    seeds: anaerobic_codh
    primary_pfam: PF03598   # CODH_C, anaerobic Ni-Fe CODH
    supporting_pfam:
      - PF02906             # Fe_hyd_lg
    hmm_profiles:
      - data/profiles/cooS.hmm
    seed_uniprot: [P31896]
window_genes: 10
known_organisms_table: data/references/known_co_organisms.csv
expected_sensory_chemistry:
  - heme
  - iron_sulfur_4Fe4S
seed: 20260429
```

`cn.yaml`:

```yaml
analyte: CN
display_name: "hydrogen cyanide and cyanide-tolerant respiration"
anchor_seeds:
  cyd_oxidase: data/seeds/cn_cyd_seeds.faa
  hcn_synthase: data/seeds/cn_hcn_synthase_seeds.faa
anchor_subfamilies:
  - name: cyd_oxidase
    seeds: cyd_oxidase
    primary_pfam: PF01654   # cydA
    supporting_pfam:
      - PF02322             # cydB / Cyt_bd_oxidase
    hmm_profiles:
      - data/profiles/cydA.hmm
      - data/profiles/cydB.hmm
    seed_uniprot: [P0ABJ9, P0ABK2]
    accessory_rules:
      cydX_like:
        max_aa_length: 60
        co_localization_nt: 200
        same_strand: true
        adjacent_to: [cydA, cydB]
  - name: hcn_synthase
    seeds: hcn_synthase
    primary_pfam: PF00890   # FAD_binding (HcnA family)
    supporting_pfam: []
    hmm_profiles:
      - data/profiles/hcnA.hmm
      - data/profiles/hcnB.hmm
      - data/profiles/hcnC.hmm
    seed_uniprot: []   # populate from P. aeruginosa hcnABC
anchor_family_groups:
  cyd_oxidase_group:
    members: [cydA, cydB, appC, appB, cydX]
    notes: "Cytochrome bd-I and bd-II share core regulatory architecture."
window_genes: 10
known_organisms_table: data/references/known_cn_organisms.csv
expected_sensory_chemistry:
  - iron_sulfur_4Fe4S
  - cysteine_metal
  - heme
seed: 20260429
```

The HMM profile files at `data/profiles/*.hmm` are built by the
`gasregnet build-profiles` subcommand introduced in section 9, which
downloads canonical seed alignments from Pfam and builds per-anchor
profiles with `hmmbuild`.

### 6.4 `configs/regulator_families.yaml` (additions and reorderings)

Additions to V1:

```yaml
regulator_families:
  # V1 entries retained.
  - family: NorR
    class: sigma54_activator
    pfam_required: ["PF00158"]   # Sigma-54 interaction domain
    pfam_optional: ["PF02954"]   # HTH_8 (NorR HTH)
  - family: CRP_FNR_extended
    class: one_component
    pfam_required: ["PF00027"]
    pfam_optional: ["PF00325"]
    notes: "Includes CooA, FNR, ANR. Heme vs Fe-S resolved by paired_evidence."
  - family: CooA_like
    class: one_component
    pfam_required: ["PF00027", "PF00325"]
    pfam_optional: []
    sensor_role_default: sensor
    notes: "CooA-class regulators with N-terminal heme. The protein is both regulator and sensor; sensor_role assigned to the same row."
  - family: ArsR_SmtB
    class: one_component
    pfam_required: ["PF01022"]
    pfam_optional: []
    sensor_role_default: sensor
    notes: "Metalloregulator that binds metal ions through cysteine clusters."

# Resolution order: a gene matches all entries it satisfies. Final
# regulator_class is the family-class with the highest accumulated weight,
# where:
#   pfam_required all matched: +2 weight
#   pfam_optional any matched: +0.5 each
#   ties broken by family priority order (CooA_like > NorR > FNR > others).
```

The classifier in `gasregnet.annotation.regulators` is rewritten to
implement the accumulating-evidence rule. See section 7.2.

### 6.5 `configs/scoring_v2.yaml`

Replaces `configs/scoring.yaml`. The weights are the V1 starting point.
After the Bayesian fit in T-V12, the headline run uses the cross-validated
weights stored at `configs/scoring_fitted.yaml`, and the manifest records
which weight set was used.

```yaml
locus_score_weights:
  anchor_marker: 3.0
  accessory_marker: 1.0
  operon_integrity: 1.5
  homology_confidence: 1.5         # raised; profile-based homology is now real
  taxonomic_context: 1.0           # raised; ecology table is now wired
  neighborhood_completeness: 0.5
  conservation_across_taxa: 1.0    # NEW

candidate_score_weights:
  locus: 1.0
  regulator_domain: 1.5
  sensor_domain: 2.0               # renamed from sensory_domain
  proximity: 1.0
  orientation: 0.25
  archetype_conservation: 1.5
  enrichment: 2.0
  taxonomic_breadth: 0.5
  structural_plausibility: 1.5
  embedding_similarity: 1.5        # NEW (P1)
  foldseek_similarity: 1.0         # NEW (P1)
  operator_motif: 2.0              # NEW
  paired_sensor_evidence: 2.0      # NEW

confidence_thresholds:
  high: 8.0    # raised
  medium: 5.0
  low: 2.0

enrichment:
  test: cmh                          # Cochran-Mantel-Haenszel
  fallback: fisher
  multiple_comparison: benjamini_hochberg
  alpha: 0.05
  case_control_ratio: [1, 3]
  permutations: 10000
  stratify_by: genus                 # NEW
  deduplication_robustness:
    - none
    - one_per_genus
    - one_per_family

windows:
  strict: 5
  standard: 10
  extended: 20
  default: standard

robustness:
  windows_to_test: [5, 10, 20]
  weight_perturbation_pct: 20
  bootstrap_iterations: 1000          # NEW

bayesian_fit:
  enabled: false                      # set to true when the headline run is fit
  model: pymc_logistic_regression
  n_samples: 2000
  n_tune: 1000
  cv_folds: 5
  cv_strategy: leave_one_genus_out
```

## 7. New module specifications

Each subsection lists module purpose, signatures, behavior, and tests.
Function bodies are the agent's responsibility. Every new module gets a
Pandera schema if it produces a cross-module table.

### 7.1 `gasregnet/search/anchors.py` (new module)

Replaces term-scan as the anchor-detection primary path.

```python
def detect_anchors_profile(
    catalog_db: Path,
    *,
    analyte: AnalyteConfig,
    profile_dir: Path,
    diamond_db: Path | None,
    seed_faa: Path,
    bitscore_threshold: float | None = None,
    e_value_threshold: float = 1e-20,
) -> pl.DataFrame:
    """Run pyhmmer profile search and DIAMOND back-search, return AnchorHits."""
```

Behavior. For each anchor subfamily, run pyhmmer's HMMER backend on the
catalog's translated CDS FASTA. Apply the curator-set bitscore threshold
or the per-profile gathering threshold. For each hit, run DIAMOND from
the hit sequence back against `seed_faa` and keep only hits whose best
DIAMOND match is to the same subfamily's seed set. Emit `AnchorHits`
schema with `evidence_type='profile_match'`. Term-scan results from V1
remain available as `evidence_type='term_match'` and are subordinated
through a `confidence` column that reads `profile_match` higher than
`term_match`.

Tests. Synthetic catalog with one true coxL, one cooS, one off-target
xanthine dehydrogenase. The function returns the coxL under `aerobic_codh`
and the cooS under `anaerobic_codh`, and the xanthine dehydrogenase is
absent at default thresholds.

### 7.2 `gasregnet/annotation/roles.py` (new module)

Three-role assignment supersedes V1's regulator-only logic.

```python
def assign_sensor_roles(
    genes: pl.DataFrame,
    *,
    regulator_families: list[RegulatorFamilyEntry],
    sensory_domain_catalog: list[SensoryDomainEntry],
    paired_evidence_rules: list[PairedEvidenceRule],
) -> pl.DataFrame:
    """Add sensor_role and refined regulator_class columns to genes."""
```

Behavior. For each gene, compute three boolean evidence vectors. The
anchor vector is true when the gene is the anchor of its locus. The
regulator vector is set by the accumulating-evidence rule from section
6.4. The sensor vector is set when at least one matched Pfam has
`role=sensor` in the sensory-domain catalog after paired-evidence
rescoring. The `sensor_role` enum is computed deterministically:

```text
anchor:    is_anchor == True
regulator: regulator_evidence > 0 and is_anchor == False
sensor:    sensor_evidence > 0 and regulator_evidence == 0 and is_anchor == False
both:      sensor_evidence > 0 and regulator_evidence > 0 and is_anchor == False
accessory: in same locus, not anchor/regulator/sensor, has any annotated function
none:      no annotation
```

`sensor_role='both'` covers CooA, RcoM, FNR, and ArsR-class regulators
where the regulator and the sensor are the same protein. The candidate
scoring picks up `sensor_role='both'` and `sensor_role='sensor'` for
sensor scoring; regulator scoring counts both `regulator` and `both`.

A new `sensor_regulator_pairs` table records two-component pairs in
the same locus, linking the HK and the RR by intergenic distance and
strand. See section 8.

Tests. CooA classified `both`; RcoM classified `both`; FixL classified
`sensor`; FixJ in the same locus classified `regulator`; the pair
appears in `sensor_regulator_pairs`.

### 7.3 `gasregnet/operator/motifs.py` (new module)

```python
def extract_upstream_regions(
    loci: pl.DataFrame,
    catalog_db: Path,
    *,
    upstream_nt: int = 200,
    downstream_nt: int = 50,
) -> pl.DataFrame:
    """Pull anchor upstream sequences from contig FASTA."""

def discover_motifs_streme(
    upstream_fasta: Path,
    *,
    out_dir: Path,
    n_motifs: int = 5,
    seed: int,
) -> pl.DataFrame:
    """Run STREME and parse motif PWMs."""

def match_collectf_regprecise(
    motifs: pl.DataFrame,
    references: Path,
) -> pl.DataFrame:
    """Match discovered motifs to known TF binding-site PWMs."""

def score_operator_support(
    candidates: pl.DataFrame,
    matches: pl.DataFrame,
) -> pl.DataFrame:
    """Add operator_motif_score per candidate."""
```

Behavior. Extract upstream regions for high-confidence loci per analyte
into one FASTA. Run STREME from the MEME suite with seed pinned for
reproducibility. Tomtom or direct PWM cosine similarity matches against
CollecTF and RegPrecise PWMs that are curated for one-component
regulators. Per-candidate score is the maximum match score over all
discovered motifs whose PWM matches the candidate's regulator family at
q < 0.1. The score contributes to `candidate_score` through the new
weight `operator_motif`.

Tests. Synthetic upstream sequence with embedded FNR consensus and a
candidate annotated as FNR class produces a non-zero motif score; the
same upstream with the consensus removed produces zero.

### 7.4 `gasregnet/scoring/conservation.py` (new module)

Replaces V1's `archetype_conservation_score = 0.0` placeholder.

```python
def compute_conservation_scores(
    candidates: pl.DataFrame,
    archetypes: pl.DataFrame,
    loci: pl.DataFrame,
    *,
    min_loci_per_archetype: int = 3,
) -> pl.DataFrame:
    """Compute archetype_conservation_score and taxonomic_breadth_score per candidate."""
```

Behavior. For each candidate, look up its archetype membership. Score
zero if the archetype has fewer than `min_loci_per_archetype`. Otherwise
score the geometric mean of three quantities clipped to `[0, 1]`:

- `genus_breadth = (n_distinct_genera_in_archetype) / (n_loci_in_archetype)`
- `position_conservation = (count_of_loci_in_archetype_with_same_regulator_class_at_same_relative_index) / (n_loci_in_archetype)`
- `family_breadth = (n_distinct_families_in_archetype) / (n_loci_in_archetype)`

This penalizes archetypes dominated by one clonal cluster while rewarding
deep cross-genus and cross-family conservation of the regulator's
position and class. Additionally compute `taxonomic_breadth_score` per
candidate as `n_distinct_phyla_in_archetype / 30` clipped to 1.

Tests. An archetype with twelve loci spanning seven genera scores
substantially higher than one with twelve loci all from the same genus.

### 7.5 `gasregnet/scoring/enrichment.py` (rewrite)

Replaces V1 Fisher's-only enrichment.

```python
def run_stratified_enrichment(
    case_genes: pl.DataFrame,
    control_genes: pl.DataFrame,
    *,
    analyte: str,
    case_definition: str,
    control_definition: str,
    stratum_column: str = "genus",
    test: str = "cmh",
    alpha: float = 0.05,
) -> pl.DataFrame:
    """Cochran-Mantel-Haenszel stratified by genus, BH corrected."""
```

Behavior. For each feature, build a stratum-by-2x2 contingency stack:
one 2x2 table per genus. Compute the CMH common odds ratio with
Mantel-Haenszel weights through statsmodels'
`StratifiedTable.test_null_odds()`. Report stratum-specific odds ratios
in a side table. BH correct across features per analyte.

Add a deduplication-robustness panel that re-runs the same enrichment
under three policies: no deduplication, one representative per genus
(highest locus_score), one representative per family (highest
locus_score). Persist the three sets of q-values in a long-format
`enrichment_robustness` table. The candidate `enrichment_score` uses the
strictest policy as default; the manifest records which.

Tests. Twenty Pseudomonas case loci and twenty matched controls with
five genera contributing 4 case+4 control each: under V1 Fisher's, FNR/CRP
appears enriched at small p-value; under CMH, the enrichment is correctly
attributed only to genera where it actually concentrates.

### 7.6 `gasregnet/structure/align.py` (new module, replaces zip-by-order)

```python
def tmalign(
    model_pdb: Path,
    homolog_pdb: Path,
    *,
    chain_id_model: str | None = None,
    chain_id_homolog: str | None = None,
    binary: str = "USalign",
) -> StructureAlignmentResult:
    """Run US-align (preferred) or TM-align and parse the alignment."""

@dataclass
class StructureAlignmentResult:
    tm_score_query: float
    tm_score_target: float
    rmsd: float
    aligned_length: int
    seq_identity_aligned: float
    residue_mapping: pl.DataFrame   # model_residue_number, homolog_residue_number
    rotation_matrix: np.ndarray     # 3x3
    translation_vector: np.ndarray  # 3
```

Behavior. Subprocess to USalign if present; fall back to TMalign.
Parse the alignment output for the residue-residue mapping and the
transformation matrix. Where both binaries are absent and the agent is
running in a constrained environment, fall back to a Biopython-based
sequence-anchored superposition: align sequences with MMseqs2 or MAFFT,
extract paired CA atoms, run `Bio.PDB.Superimposer`, and emit a
`StructureAlignmentResult` with `tm_score_query=NaN` and a manifest
warning. Never zip CA atoms by raw position order.

Tests. Two homologous PDB chains of different lengths produce a residue
mapping where conservation of a known catalytic site is correctly
identified.

### 7.7 `gasregnet/embeddings/esm.py` (new module, P1)

```python
def embed_proteins_esm2(
    sequences: pl.DataFrame,
    *,
    model_name: str = "esm2_t33_650M_UR50D",
    batch_size: int = 8,
    device: str = "cpu",
    cache_dir: Path,
) -> pl.DataFrame:
    """Compute mean-pooled per-protein embeddings."""

def nearest_known_sensor(
    candidate_embeddings: pl.DataFrame,
    benchmark_embeddings: pl.DataFrame,
    *,
    metric: str = "cosine",
) -> pl.DataFrame:
    """For each candidate, find nearest positive-control sensor and distance."""
```

Behavior. Use `fair-esm` with the t33 650M model by default. The 3B model
is gated on the GPU extra. Cache embeddings under `cache_dir` keyed by
SHA-256 of the sequence and the model name; subsequent runs are
near-instant. Output schema: `gene_accession`, `model_name`, `dim` (1280
for t33 650M), `embedding` (List[Float32]).

The `embedding_similarity_score` in candidate scoring is
`max(0, 1 - cosine_distance(candidate, nearest_positive_in_benchmark))`.
Independent of sequence-level similarity, this finds candidates that
fold into the same family.

Tests. Two known CooA paralogs from different organisms have cosine
similarity > 0.95 in ESM-2 space; a CooA paralog and a generic LacI have
similarity < 0.6. Use a tiny test fixture and `esm2_t6_8M_UR50D` for the
unit test to keep CI fast.

### 7.8 `gasregnet/structure/foldseek.py` (new module, P1)

```python
def foldseek_search_afdb(
    query_pdb: Path,
    *,
    afdb_index: Path,
    out_tsv: Path,
    e_value: float = 1e-3,
    max_results: int = 100,
) -> pl.DataFrame:
    """Run Foldseek easy-search against an AlphaFold DB index."""
```

Behavior. Subprocess to the foldseek binary. Parse the output into a
table with columns `query`, `target`, `tm_score`, `lddt`, `e_value`,
`prob`, `target_organism`. The candidate's
`foldseek_similarity_score` is the maximum prob over all hits whose
target is annotated as a known sensor in the benchmark or a structural
neighbor of one.

The AFDB index is not bundled. The CLI subcommand
`gasregnet build-foldseek-index --afdb-bacterial` builds a local
Foldseek index from the AFDB bacterial subset (or a curated subset
defined by GTDB phylum coverage).

Tests. Mock the foldseek binary in CI; integration test uses a fixture
TSV.

### 7.9 `gasregnet/scoring/bayesian.py` (new module, P1)

```python
def fit_bayesian_scoring(
    candidates: pl.DataFrame,
    benchmark_recovery: pl.DataFrame,
    *,
    feature_columns: list[str],
    cv_strategy: str = "leave_one_genus_out",
    n_samples: int = 2000,
    n_tune: 1000,
    seed: int,
) -> BayesianScoringResult:
    """Fit a Bayesian logistic regression to predict known-sensor labels."""

@dataclass
class BayesianScoringResult:
    posterior_weights: pl.DataFrame  # mean and 94% HDI per feature
    cv_metrics: pl.DataFrame         # AUC, AUPRC, recall@5, recall@10 per fold
    weights_yaml: Path
```

Behavior. Build a feature matrix where columns are the named score
components (regulator_domain, sensor_domain, proximity, conservation,
enrichment, embedding_similarity, etc.) and rows are
benchmark-positive plus matched-negative candidates. Fit a Bayesian
logistic regression in pymc with weakly informative priors centered on
zero. Run leave-one-genus-out CV. Report the posterior mean weights and
their 94% HDIs. Write the fitted weights to
`configs/scoring_fitted.yaml` and the headline run reads from there.

The fast baseline is scikit-learn's `LogisticRegressionCV` with
`StratifiedGroupKFold(group=genus)`. The pymc model is the headline. The
agent implements both and the CLI exposes a flag to choose.

Tests. On a synthetic benchmark where the true generative weights are
known, the recovered posterior means lie within their HDIs around the
true weights at the 1000-sample scale.

### 7.10 `gasregnet/db/duckdb_views.py` (new module)

DuckDB analytical kernel.

```python
def build_results_database(
    results_dir: Path,
    *,
    out_db: Path,
) -> Path:
    """Build a single .duckdb file exposing every Parquet result as a view."""
```

Behavior. Open `out_db`, register Parquet files at
`results_dir/intermediate/*.parquet` as DuckDB views, define curator-named
views that join loci, genes, candidates, archetypes, and enrichment
through the foreign keys in the schemas. Persist the view DDL into the
database. After this runs, a biologist can type:

```sql
duckdb results/headline/results.duckdb
> select c.candidate_id, c.regulator_class, c.candidate_score, l.organism, a.architecture_string
  from candidates c
  join loci l using (locus_id)
  left join archetypes a using (archetype_id)
  where c.analyte = 'CO' and c.candidate_score > 8
  order by c.candidate_score desc
  limit 20;
```

Required views:

```text
view: high_confidence_loci      -> loci where locus_confidence in ('high','medium')
view: top_candidates             -> candidates ordered by candidate_score
view: candidates_with_locus      -> candidates JOIN loci
view: candidates_full            -> candidates JOIN loci JOIN archetypes JOIN enrichment
view: archetype_summary          -> archetypes with member loci counts and dominant features
view: enrichment_significant     -> enrichment_results where q_value < 0.05
view: benchmark_recovery_full    -> benchmark JOIN candidates JOIN loci
view: sensor_regulator_pairs     -> the new pairs table from section 7.2
```

Tests. After `make repro`, the test opens the produced DuckDB file,
runs each view definition, and asserts row counts are non-zero where
data is available.

The headline run produces both `results/headline/intermediate/*.parquet`
(durable) and `results/headline/results.duckdb` (analytical kernel).
Biologists download the .duckdb file directly and explore.

### 7.11 `gasregnet/provenance/prov.py` (new module, P1)

```python
def emit_prov_for_candidate(
    candidate_id: str,
    candidates: pl.DataFrame,
    loci: pl.DataFrame,
    enrichment: pl.DataFrame,
    archetypes: pl.DataFrame,
    *,
    manifest: dict[str, object],
    out_path: Path,
) -> Path:
    """Emit PROV-O JSON-LD for one candidate's full derivation."""
```

Behavior. For every candidate, write a JSON-LD document describing the
candidate as a `prov:Entity` and the score components as
`prov:wasDerivedFrom` other entities (the locus, the enrichment row, the
archetype membership, the embedding match). Each derivation is annotated
with the tool version, config hash, and run timestamp from the manifest.
Reviewers reconstruct exactly how each candidate was derived from a
single JSON-LD file.

Tests. Round-trip: emit the PROV doc, parse it with `prov`, and confirm
the derivation graph contains every score-component edge.

## 8. Updated data contracts

V2 extends V1 schemas with new optional columns and adds new tables.

### 8.1 Extensions to existing schemas

`GenesSchema`: add `sensor_role` (Utf8 in {anchor, regulator, sensor,
both, accessory, none}); add `paired_evidence_applied` (List[Utf8]);
add `embedding_id` (Utf8, nullable, foreign key to `embeddings`).

`RegulatorCandidatesSchema`: add `sensor_role` (Utf8); add
`embedding_similarity_score` (Float64); add `foldseek_similarity_score`
(Float64); add `operator_motif_score` (Float64); add
`paired_sensor_evidence_score` (Float64); add `bayesian_posterior_mean`
(Float64, nullable until the Bayesian fit completes); add
`bayesian_posterior_hdi_low` (Float64, nullable);
`bayesian_posterior_hdi_high` (Float64, nullable).

`LociSchema`: add `phylum`, `class`, `order`, `family`, `genus` (Utf8,
all populated by `annotate_taxonomy` before scoring); add
`conservation_across_taxa_score` (Float64).

### 8.2 New schemas

`AnchorHitsSchema` extension: add `evidence_type` enum
{profile_match, term_match, diamond_back_confirmed}; add `bitscore`,
`e_value`, `identity`, `coverage` as Float64.

`SensorRegulatorPairsSchema`:

```text
pair_id                    Utf8 primary key
analyte                    Utf8
locus_id                   Utf8 FK
hk_gene_accession          Utf8
rr_gene_accession          Utf8
intergenic_distance_nt     Int64
co_strand                  Boolean
hk_sensor_domains          List[Utf8]
hk_sensory_chemistries     List[Utf8]
rr_dna_binding_domains     List[Utf8]
pair_score                 Float64
```

`OperatorMotifsSchema`:

```text
motif_id                   Utf8 primary key
analyte                    Utf8
discovery_method           Utf8 (streme, rsat, regprecise_match)
n_sites_observed           Int32
match_collectf_id          Utf8 nullable
match_regprecise_id        Utf8 nullable
match_q_value              Float64 nullable
pwm                        Utf8 (MEME-formatted PWM as text blob)
```

`EmbeddingsSchema`:

```text
embedding_id               Utf8 primary key
gene_accession             Utf8
model_name                 Utf8
dim                        Int32
sequence_sha256            Utf8
embedding                  List[Float32]   # length == dim
```

`StructureAlignmentResultsSchema`:

```text
alignment_id               Utf8 primary key
candidate_id               Utf8 FK
homolog_pdb_id             Utf8
tm_score_query             Float64
tm_score_target            Float64
rmsd                       Float64
aligned_length             Int32
seq_identity_aligned       Float64
n_conserved_residues       Int32
binding_pocket_residues    List[Int32]
```

`EnrichmentRobustnessSchema`:

```text
analyte                    Utf8
deduplication_policy       Utf8 (none, one_per_genus, one_per_family)
feature_type               Utf8
feature_name               Utf8
n_case                     Int64
n_control                  Int64
odds_ratio                 Float64
p_value                    Float64
q_value                    Float64
concordance_with_default   Float64    # Spearman with the default policy
```

## 9. Updated CLI

Subcommand additions. Every new subcommand follows V1 conventions
(`--config`, `--out`, `--seed`, `--verbose`, manifest emission,
config.resolved.yaml).

```bash
gasregnet check-tools                         # probe external binaries, write tools_resolved.yaml
gasregnet build-profiles                      # download Pfam alignments, hmmbuild for anchor families
gasregnet detect-anchors-profile              # profile-driven anchor detection (replaces smoke as default)
gasregnet assign-roles                        # three-role gene annotation
gasregnet extract-upstream                    # pull anchor upstream sequences
gasregnet discover-motifs                     # STREME on high-confidence loci
gasregnet match-motifs                        # against CollecTF and RegPrecise
gasregnet score-conservation                  # compute archetype conservation
gasregnet enrich-stratified                   # CMH + dedup robustness
gasregnet align-structure                     # USalign / TMalign wrapper
gasregnet embed-proteins                      # ESM-2 embeddings
gasregnet foldseek-search                     # against AFDB index
gasregnet build-foldseek-index                # one-time AFDB ingest
gasregnet fit-bayesian                        # Bayesian scoring fit
gasregnet build-views                         # DuckDB analytical kernel
gasregnet emit-provenance                     # PROV-O for top candidates
```

The headline pipeline is:

```bash
make sync
make assets
make profiles      # NEW: build HMM profiles
make datasets
make index-datasets
make detect-anchors-profile     # replaces detect-anchors term-scan
make assign-roles
make extract-neighborhoods
make extract-upstream
make discover-motifs
make match-motifs
make score
make score-conservation
make enrich-stratified
make embeddings
make foldseek
make fit-bayesian
make report
make build-views
make emit-provenance
make repro
```

Each Make target invokes the corresponding Snakemake rule from
`workflows/full_v2.smk`. A new top-level workflow file consolidates the
V2 DAG.

## 10. Task graph for V2

V2 tasks continue the V1 numbering. All V1 tasks T0 through T18 stand.
V2 introduces T-V1 through T-V14.

Every V2 task has a definition of done. Where a task corrects a V1
defect, the defect number from section 2 is in parentheses.

### T-V1. Replace benchmark with regulator benchmark v2

Defect 1. Owner: 1 agent. Depends on T-V0 (none beyond V1).

Output: `configs/benchmarks/regulators_v2.csv`,
`data/benchmarks/regulators_v2.csv`, `data/benchmarks/anchors_v1.csv` (V1
file renamed in place). Update `configs/benchmarks.yaml` to reference
both. Update `BenchmarkSchema` if necessary (it should not be).

DOD: every entry in the table from section 6.1 is present with UniProt
accession and first-publication citation. `verify_pmid` column flags any
PMID the agent could not directly confirm. Schema validates. Unit test
asserts at least 10 entries per analyte plus 10 negative controls.

### T-V2. Correct sensory_domains.yaml

Defect 3. Owner: 1 agent. Depends on T-V1.

Output: `configs/sensory_domains.yaml` rewritten per section 6.2,
including the `paired_evidence` block. V1 file moves to
`configs/legacy/v1/sensory_domains.yaml`.

DOD: schema validates with new fields (`role`, `notes`, `paired_evidence`).
Update `SensoryDomainEntry` Pydantic model. All callers of the old config
load through the new model. Coverage on `gasregnet/config.py` remains
above 95%.

### T-V3. Split CO seeds and rebuild HMM profiles

Defects 2, 13. Owner: 1 agent. Depends on T-V2.

Output: `data/seeds/co_aerobic_seeds.faa`,
`data/seeds/co_anaerobic_seeds.faa`,
`data/seeds/cn_cyd_seeds.faa`, `data/seeds/cn_hcn_synthase_seeds.faa`.
A new `scripts/build_profiles.py` that uses pyhmmer to build HMMs from
seed alignments downloaded from the Pfam/InterPro APIs. Outputs at
`data/profiles/*.hmm`.

DOD: every anchor family in `co.yaml` and `cn.yaml` has a corresponding
profile. `pyhmmer.plan7.HMMFile.read` loads each. SHA-256 of each
profile is recorded in `configs/assets.yaml`.

### T-V4. Profile-driven anchor detection module

Defect 2. Owner: 1 agent. Depends on T-V3.

Output: `gasregnet/search/anchors.py` per section 7.1. Wire
`gasregnet detect-anchors-profile` into the CLI.

DOD: synthetic catalog test in section 7.1 passes. Term-scan retained as
fallback under `--mode smoke`. The default mode is `profile`.

### T-V5. Sensor role classification

Defects 4, 5. Owner: 1 agent. Depends on T-V2.

Output: `gasregnet/annotation/roles.py` per section 7.2. New
`SensorRegulatorPairsSchema`. Update `GenesSchema` per section 8.1.

DOD: CooA classified `sensor_role='both'`; FixL classified `sensor`;
FixJ classified `regulator`; the FixL+FixJ pair appears in
`sensor_regulator_pairs`. Unit tests cover three more known pairs.

### T-V6. Operator and motif module

Defect 14. Owner: 1 agent. Depends on T-V4.

Output: `gasregnet/operator/motifs.py` per section 7.3. New
`OperatorMotifsSchema`. Bundle a small subset of CollecTF/RegPrecise
PWMs at `data/references/binding_sites/`. The `make match-motifs` target
runs end-to-end on the corpus mode fixture.

DOD: synthetic upstream test in section 7.3 passes. The candidate
scorer reads `operator_motif_score` from the new column.

### T-V7. Conservation layer

Defect 8. Owner: 1 agent. Depends on T-V5, V1's archetype clustering.

Output: `gasregnet/scoring/conservation.py` per section 7.4. Wire into
the candidate scorer.

DOD: synthetic test in section 7.4 passes. Headline run produces
non-zero `archetype_conservation_score` for at least 25 percent of
candidates.

### T-V8. Stratified enrichment

Defect 9. Owner: 1 agent. Depends on T-V5.

Output: rewrite `gasregnet/scoring/enrichment.py` per section 7.5. New
`EnrichmentRobustnessSchema`. Update existing tests to compare CMH and
Fisher results on the synthetic-stratification fixture.

DOD: synthetic clade-imbalance test in section 7.5 passes. Robustness
panel emits three q-value sets.

### T-V9. Wire operon integrity and ecology into ingest

Defect 7. Owner: 1 agent. Depends on V1's neighborhoods/operons.py.

Output: `gasregnet/io/sqlite_efi.py` and
`gasregnet/datasets/refseq.py` updated to call
`anchor_operon_integrity` and `score_taxonomic_context` before emitting
loci rows.

DOD: emitted loci have non-zero `operon_integrity_score` whenever a real
operon co-orientation exists. Coverage of the changed lines is 100
percent.

### T-V10. Real structural alignment

Defect 6. Owner: 1 agent. Depends on V1.

Output: `gasregnet/structure/align.py` per section 7.6. Mark
`gasregnet/structure/pdb.residue_mapping_by_order` as deprecated and
raise `DeprecationWarning` on call. Wire `align.tmalign` into Figure 6
generation.

DOD: integration test against two real homologous PDB chains gives a
non-degenerate alignment with TM-score above 0.5 and a residue-residue
mapping that respects the structural alignment, not the sequence index.

### T-V11. ESM-2 embedding layer

Defect 14, P1. Owner: 1 agent. Depends on T-V1.

Output: `gasregnet/embeddings/esm.py` per section 7.7. New
`EmbeddingsSchema`. New CLI subcommand `gasregnet embed-proteins`. Add
`fair-esm` and `torch` to the optional `embeddings` extra.

DOD: CooA-paralog cosine similarity test passes. Cache hit on second
run is sub-second per protein. Embedding cache stored at
`cache/embeddings/<model_name>/<sha256>.npy`.

### T-V12. Foldseek integration

P1. Owner: 1 agent. Depends on T-V10.

Output: `gasregnet/structure/foldseek.py` per section 7.8. New CLI
subcommands `gasregnet build-foldseek-index` and `foldseek-search`.

DOD: CI runs against a mock foldseek binary; integration test loads a
fixture TSV. Manifest records foldseek version.

### T-V13. Bayesian scoring

P1. Owner: 1 agent. Depends on T-V1, T-V11, T-V8.

Output: `gasregnet/scoring/bayesian.py` per section 7.9. CLI
`gasregnet fit-bayesian`. Headline `configs/scoring_fitted.yaml`
written by the fit. Manifest records which weight set the report run
used.

DOD: synthetic-recovery test passes. Real benchmark fit has at least
five-fold leave-one-genus-out CV with reported AUC and AUPRC. The
posterior weights file is loadable by V1's `ScoringConfig` model.

### T-V14. DuckDB analytical kernel

Defect: V1's analytical layer is Polars-only. Owner: 1 agent. Depends
on V1, all V2 schema additions stable.

Output: `gasregnet/db/duckdb_views.py` per section 7.10. CLI
`gasregnet build-views`. Headline run produces
`results/headline/results.duckdb` alongside the Parquet intermediates.

DOD: every view in section 7.10 returns non-zero rows on the headline
fixture. README has a "Querying results" section showing the example
SQL above.

### T-V15. Provenance emission

P1. Owner: 1 agent. Depends on T-V14.

Output: `gasregnet/provenance/prov.py` per section 7.11. CLI
`gasregnet emit-provenance`.

DOD: round-trip test in section 7.11 passes. Headline run produces
PROV-O JSON-LD for the top 30 candidates per analyte.

## 11. Definition of done at V2

The codebase is V2-complete and preprint-ready when all of the following
hold.

1. `make repro` runs end to end in under 90 minutes on a developer laptop
   using the SQLite fixture and produces the six headline figures, the
   T1-T6 tables, the partition-outcome JSON, the DuckDB analytical
   database, the Foldseek hit table for the top 30 candidates, the
   Bayesian fit weights, and the PROV-O JSON-LD for the top 30
   candidates.
2. The benchmark recovery integration test, run against `regulators_v2`,
   reports recall above 0.6 on the positives and false positive rate
   below 0.1 on the negatives.
3. The partition claim integration test reports the chemistry partition
   outcome and writes `partition_outcome.json` with q-value, contingency
   table, and the chemistry-stratified contingency table.
4. The robustness panel reports at least 0.7 Spearman concordance of the
   top-100 candidate ranking across the three deduplication policies and
   across windows {5, 10, 20}.
5. Every score component column is populated for every candidate; no
   candidate has all-zero score components.
6. `mypy --strict` is clean. `ruff check` is clean. Coverage on
   `gasregnet/` is 80 percent or higher.
7. The DuckDB analytical kernel exposes every view in section 7.10 and a
   biologist can run the example SQL query from section 7.10 against the
   headline run.
8. The PROV-O documents for the top candidates can be parsed and the
   derivation graph contains every score-component edge.
9. `manifest.json` records the package version, every config hash, every
   input data hash, the seed, the wall clock, the version of every
   external binary including HMMER, DIAMOND, MMseqs2, MEME, Foldseek,
   USalign, and the Bayesian-fit weights file path and SHA-256.
10. The repository has a tagged release at `v0.2.0` with a Zenodo DOI
    and a CITATION.cff that includes the V2 contributors and binaries.

When all ten conditions are met, the codebase is the preprint's Methods
section made executable for the V2 manuscript.

## 12. Anti-patterns (additions to V1)

Do not run anchor detection by substring matching except as a permissive
recall-augmenting layer with `evidence_type='term_match'`.

Do not score sensors and regulators in one component. They are different
biology and the partition claim depends on keeping them separate.

Do not list HAMP, GGDEF, EAL, or HD-GYP as sensory domains. They are
transducers and effectors. Do not let them contribute to
`sensor_domain_score`.

Do not zip CA atoms by sequence order and call it structural alignment.
Use TM-align or US-align.

Do not run Fisher's exact on phylogenetically clumped loci without
stratification. CMH is the default. Fisher's is the fallback.

Do not commit a benchmark of the genes the pipeline already searches
for. The benchmark must be the regulators the pipeline is supposed to
find.

Do not compute distances in gene-index space and call them nucleotide
distances. Maintain `relative_index` and `distance_nt` separately, and
populate `distance_nt` from contig coordinates only.

Do not run the embedding model in inference loops without caching.
ESM-2 inference is the slowest step; cache by sequence SHA-256 and model
name.

Do not write SQL views that hide their join conditions. Every view in
the analytical kernel has a comment block explaining its joins and its
intended biologist use case.

Do not skip the PROV-O emission for the top candidates. Reviewers can
and will reconstruct derivations from it.

## 13. Parallelization plan

Five-agent slice for V2 after the configs land (T-V1, T-V2, T-V3
sequential; one agent each):

- Agent A: T-V4 (anchor detection) then T-V8 (stratified enrichment)
- Agent B: T-V5 (roles) then T-V7 (conservation)
- Agent C: T-V6 (motifs) then T-V13 (Bayesian fit)
- Agent D: T-V9 (operon and ecology wiring) then T-V10 (structural alignment)
- Agent E: T-V11 (ESM-2) then T-V12 (Foldseek) then T-V14 (DuckDB views)

T-V15 (provenance) lands last because it depends on the analytical
kernel from T-V14.

T-V1, T-V2, T-V3 are gating. They define the contracts and assets every
downstream agent reads. Land them strictly in order, with one agent at a
time.

## 14. Final note to agents

V1 was a tool. V2 is a discovery. The biology in this document is more
important than the engineering. When in doubt, prefer biological
correctness to architectural elegance. If a Pandera schema makes a
biological refinement awkward, change the schema. If a config field name
reads cleanly in code but reads wrong to a biochemist, rename the field.

Two principles to internalize.

First. Sensors, regulators, and effectors are three different roles. The
current pipeline conflates them. Disentangling them is the load-bearing
correction for the partition claim that becomes the headline preprint
result.

Second. Co-occurrence is not regulation. The motif and operator module
in section 7.3 is what converts the candidate list from "regulators near
gas-anchor genes" into "regulators that plausibly bind operators of
gas-anchor genes". That is the difference between a tool paper and a
finding paper.

V2 is what the preprint ships on. The preprint's Discussion frames V1 as
the discovery framework and V2 as the rigor pass. The reader does not
need to know that. The reader needs to read a biologically credible
manuscript built on a reproducible pipeline whose code does what the
methods say it does. V2 makes that true.