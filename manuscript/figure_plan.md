# scRFU manuscript figure plan

## Frozen scope

This plan uses only frozen scRFU evidence. It defines a first-pass visual story;
it does not authorize new analysis or modification of validated source tables.

> scRFU provides a scalable, scverse-native functional-unit representation
> that maps heterogeneous T-cell receptor sequences into a transferable
> repertoire feature space while preserving single-cell context.

Panel sources and caveats are specified in
[`source_table_to_panel_map.tsv`](source_table_to_panel_map.tsv). Portable source
labels resolve against the sealed external evidence roots in
`docs/evidence_index.json`.

## Editorial compression

The initial plan contained 28 main panels. The compressed plan contains **22**:

| Figure | Initial | Final | Decision |
|---|---:|---:|---|
| Figure 1 | 5 | 5 | Retain five compact architecture panels |
| Figure 2 | 6 | 4 | Move detailed chunking and storage overhead to Extended Data |
| Figure 3 | 5 | 5 | Retain the complete representation argument |
| Figure 4 | 6 | 4 | Merge cohort design into captions and combine held-out coverage/stability |
| Figure 5 | 6 | 4 | Merge context schematic with phenotype panel and match coverage with coherence |

### Panel disposition

| Original panel | Classification | Final disposition |
|---|---|---|
| 1A concept | ESSENTIAL MAIN | Figure 1A |
| 1B native storage | ESSENTIAL MAIN | Split into chain alignment (1B) and storage (1C) |
| 1C workflow/API boundary | ESSENTIAL MAIN | Integrated into 1A–1C rather than retained as a separate software-manual panel |
| 1D public parity | ESSENTIAL MAIN | Figure 1E |
| 1E cell-summary policy | SUPPORTING MAIN | Figure 1D |
| 2A official parity | ESSENTIAL MAIN | Figure 2A |
| 2B full Wells | ESSENTIAL MAIN | Figure 2B |
| 2C detailed chunk/worker scaling | MOVE TO EXTENDED DATA | ED2; main invariance badge retained in 2D |
| 2D restart/cache | SUPPORTING MAIN | Figure 2D |
| 2E native 25k/100k/250k | ESSENTIAL MAIN | Figure 2C |
| 2F serialized storage cost | MOVE TO EXTENDED DATA | ED3 |
| 3A compression | ESSENTIAL MAIN | Figure 3A |
| 3B sequences per RFU | SUPPORTING MAIN | Figure 3B |
| 3C dimension/sparsity | ESSENTIAL MAIN | Figure 3C |
| 3D coverage | ESSENTIAL MAIN | Figure 3D |
| 3E feature sharing | ESSENTIAL MAIN | Figure 3E |
| 4A cohort-role schematic | REDUNDANT | Reduced to caption/header context |
| 4B retrieval | ESSENTIAL MAIN | Figure 4B |
| 4C within/between | ESSENTIAL MAIN | Figure 4A |
| 4D GSE190905 downsampling | SUPPORTING MAIN | Figure 4C |
| 4E held-out coverage | ESSENTIAL MAIN | Merged into Figure 4D |
| 4F held-out stability | ESSENTIAL MAIN | Merged into Figure 4D |
| 5A context schematic | REDUNDANT AS STANDALONE | Inset in Figure 5A |
| 5B phenotype coupling | ESSENTIAL MAIN | Figure 5A |
| 5C phenotype stability | SUPPORTING MAIN | Figure 5B |
| 5D VDJdb match coverage | SUPPORTING MAIN | Inset in Figure 5D |
| 5E VDJdb null | ESSENTIAL MAIN | Figure 5D |
| 5F native VDJdb linkage | ESSENTIAL MAIN | Figure 5C |

No evidence is discarded. Secondary implementation and sensitivity evidence is
assigned to Extended Data.

## Figure 1 — scverse-native RFU representation

**Single message.** scRFU maps AIRR receptor chains into a frozen RFU feature
space while preserving every chain, observation and portable scverse context.

| Panel | Content | Class |
|---|---|---|
| 1A | Heterogeneous receptors → exact CDR3 queries → frozen RFU reference → reusable RFU features | ESSENTIAL MAIN |
| 1B | Chain-aligned AIRR-to-scRFU records, including explicit noneligible TRA/BCR chains and multiple TRBs | ESSENTIAL MAIN |
| 1C | AnnData/MuData storage: aligned `obsm`, portable `uns` provenance and optional `obs` summary | ESSENTIAL MAIN |
| 1D | Ambiguity-aware cell-summary decision for multi-TRB cells | SUPPORTING MAIN |
| 1E | Public `wu2020_3k` confirmation: 7,544 chains, 2,931 eligible TRBs and zero native/table mismatches | ESSENTIAL MAIN |

**Move to Extended Data.** Full serialization, subsetting and concatenation test
matrix.

**Reviewer objection.** Is scRFU a native chain-aware scverse method or merely a
table wrapper that collapses receptor structure?

**Evidence support.** Yes. Architecture behavior is covered by exact tests and
the public Scirpy smoke provides empirical native/table parity and H5MU reload.

## Figure 2 — Fidelity, deterministic scale and restartability

**Single message.** scRFU reproduces the canonical RFU result exactly and
scales deterministically from native bounded objects to the full Wells atlas.

| Panel | Content | Class |
|---|---|---|
| 2A | Official assignment/reconstruction parity matrix | ESSENTIAL MAIN |
| 2B | Full Wells funnel, deduplication, total runtime and peak RSS | ESSENTIAL MAIN |
| 2C | Native 25k/100k/250k runtime and memory with zero mismatch strip | ESSENTIAL MAIN |
| 2D | Fresh versus cached/resumed runtime with serial/parallel/chunk/order invariance badges | SUPPORTING MAIN |

**Move to Extended Data.** Detailed chunk-size/worker timing (ED2), native
serialized size and write/reload cost (ED3), and adversarial parity details
(ED1).

**Reviewer objection.** Does portability change the scientific output, and is
the method practical and recoverable at atlas scale?

**Evidence support.** Yes. Canonical parity, full Wells, native scaling, exact
configuration invariance and cache reuse are frozen. Timing remains
host-specific.

## Figure 3 — Compact, reusable frozen feature space

**Single message.** Frozen RFUs transform highly sparse exact receptor identity
into a more compact feature vocabulary that is reused across independent
datasets.

| Panel | Content | Class |
|---|---|---|
| 3A | Exact CDR3/clonotype versus nearest/threshold RFU feature counts | ESSENTIAL MAIN |
| 3B | Sequences-per-RFU ECDF and singleton/median annotations | SUPPORTING MAIN |
| 3C | Feature dimensionality versus sample-feature sparsity | ESSENTIAL MAIN |
| 3D | Frozen-threshold coverage and RFU richness across four datasets | ESSENTIAL MAIN |
| 3E | Cross-dataset exact-CDR3 versus RFU sharing | ESSENTIAL MAIN |

The Wells/GSE157007 annotation must state: **4,898 held-out RFUs are
represented in Wells, whereas 3,187 exact CDR3 sequences are shared.** It must
not describe those cohorts or receptors as biologically equivalent.

**Move to Extended Data.** Full prevalence distributions and score quantiles.

**Reviewer objection.** Is RFU just a within-dataset relabeling, or does it
provide a reusable feature vocabulary with measurable compression?

**Evidence support.** Yes as a representation-property claim. Compression is
not itself evidence of improved biological signal, and threshold failure is not
calibrated OOD probability.

## Figure 4 — Transfer, robustness and complementary comparators

**Single message.** RFUs transfer without refitting, retain repeated-donor
structure and show robustness complementary to exact CDR3, genuine Scirpy
clonotypes and conventional summaries.

| Panel | Content | Class |
|---|---|---|
| 4A | GSE190905 within- versus between-donor cosine by representation | ESSENTIAL MAIN |
| 4B | Leave-one-timepoint-out top-1, top-3 and MRR by representation | ESSENTIAL MAIN |
| 4C | GSE190905 50%/75% subsampling stability by representation | SUPPORTING MAIN |
| 4D | Preregistered GSE157007 frozen coverage/richness with held-out subsampling inset | ESSENTIAL MAIN |

Comparator order is fixed across all Figure 4 panels: RFU, exact CDR3, Scirpy
clonotype, TRBV+TRBJ, CDR3 length, diversity. RFU uses a distinct but
non-dominant color; all other methods use equally legible neutral colors.

**Move to Extended Data.** Jaccard, individual ranks, complete seed-level
observations and the full comparator matrix (ED5).

**Reviewer objection.** Does the representation transfer fairly, and are
comparator outcomes shown even when RFU is tied or beaten?

**Evidence support.** Yes for independent/held-out transfer and bounded
two-visit representation. Exact CDR3 slightly exceeds RFU MRR; V/J has the
highest MRR among the count representations; coarse summaries are most stable
partly because they discard detail.

## Figure 5 — Single-cell context and external annotation coherence

**Single message.** RFU annotations retain cell-phenotype context and show
group-level coherence of external antigen annotations without claiming antigen
specificity.

| Panel | Content | Class |
|---|---|---|
| 5A | Full Wells top-30 abundant RFU × cell-type heatmap with chain-to-phenotype context inset | ESSENTIAL MAIN |
| 5B | Phenotype-coupling cosine and dominant-phenotype agreement under cell subsampling | SUPPORTING MAIN |
| 5C | Native `wu2020_3k` AIRR chain → RFU identity → CDR3/CDR3+V evidence linkage | ESSENTIAL MAIN |
| 5D | Wells observed-versus-null same-antigen pair fraction with exact-match coverage inset | ESSENTIAL MAIN |

Top phenotype RFUs are selected solely by total abundance, with ties resolved
by RFU label. They must never be selected by phenotype specificity.

**Move to Extended Data.** All 24 matching/assignment/ambiguity combinations,
grouping baselines, evidence scores and four null models (ED6).

**Reviewer objection.** Does the representation preserve single-cell
interpretability, and is there external evidence of within-group annotation
structure?

**Evidence support.** Yes descriptively. Phenotype analyses have no cell-level
inferential p-values. VDJdb supports external annotation coherence only, not
antigen specificity, epitope prediction or antigen-recognition prediction.

## Extended Data prototypes

- **ED1:** adversarial canonical parity and reconstruction details.
- **ED2:** chunking, parallelism, order invariance and cache/restart behavior.
- **ED3:** native storage/RSS, H5AD/H5MU round-trip, subsetting and guarded
  concatenation.
- **ED4:** threshold/reference-coverage sensitivity; 0.6 remains frozen.
- **ED5:** full comparator and downsampling matrices.
- **ED6:** VDJdb match tiers, ambiguity, grouping baselines and null sensitivity.
- **ED7:** scverse dependency, CI, package-install and interoperability matrix.

Only ED3, ED5 and ED6 have immediately available plot-ready local tables for
all quantitative components. The other technical figures may use audited ledger
values or remain compact pass/fail matrices during first-pass prototyping.

## Unsupported panels removed

- New UMAP colored by RFU: no separately frozen plot-ready source table.
- Deep longitudinal persistence/expansion/contraction: GSE190905 has two visits
  and GSE345124 RFU analysis was not run.
- RFU-level antigen specificity or prediction: unsupported by VDJdb evidence.
- Calibrated out-of-distribution probability: not tested.
- Universal superiority over Scirpy clonotypes or conventional summaries:
  contradicted by the complementary comparator results.

## Experiment requirement

No new experiment is required for first-pass figure design or for the stated
central claim. A full GSE345124 analysis would be required only if the scope is
expanded to deep multi-visit longitudinal dynamics.
