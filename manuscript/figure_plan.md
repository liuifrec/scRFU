# scRFU manuscript figure architecture

## Scope and central claim

This plan uses only frozen scRFU evidence. It defines figure logic and source
relationships; it does not authorize new analysis, source-table modification,
or polished figure generation.

> scRFU provides a scalable, scverse-native functional-unit representation
> that maps heterogeneous T-cell receptor sequences into a transferable
> repertoire feature space while preserving single-cell context.

The panel-level contract is
[`source_table_to_panel_map.tsv`](source_table_to_panel_map.tsv). External paths
in that map are portable labels relative to the sealed evidence roots indexed
in `docs/evidence_index.json`; they are not developer-machine paths.

## Interpretation boundaries

- RFU is complementary to exact sequence identity, Scirpy clonotypes, V/J
  usage, CDR3 length, and diversity summaries. The figures must show where
  comparators equal or exceed RFU results.
- Cross-dataset RFU reuse demonstrates a shared frozen feature vocabulary, not
  biological equivalence between cohorts or receptors.
- VDJdb results concern coherence of external antigen annotations among
  distinct matched sequences. They do not establish antigen specificity for
  an RFU.
- Threshold qualification is a frozen assignment policy. Threshold failure is
  not a calibrated out-of-distribution probability.
- GSE190905 is a six-donor, two-visit repeated-measures demonstration. It does
  not support claims about deep longitudinal trajectories.
- Performance and memory measurements describe the audited Linux development
  host. Python-only CI portability does not imply official RFU/R execution was
  validated on macOS or Windows.

## Figure 1 — Concept and scverse-native architecture

**Single message.** scRFU adds a frozen, chain-aligned RFU representation to
standard AIRR AnnData/MuData objects while preserving exact receptor and
single-cell identities.

| Panel | Content | Evidence role | Priority |
|---|---|---|---|
| 1A | Heterogeneous CDR3 sequences mapped to a shared frozen RFU feature vocabulary while exact identities remain available | Conceptual orientation | Indispensable |
| 1B | Positional alignment of `obsm["airr"]` and `obsm["scrfu"]`, portable `uns["scrfu"]` provenance, and optional `obs` summary | Storage architecture | Indispensable |
| 1C | Scirpy AIRR object → optional chain indexing → `scrfu.tl.assign_rfu` → explicit summary → Scanpy/Scirpy downstream use | Ecosystem and responsibility boundary | Indispensable |
| 1D | Public `wu2020_3k` native-versus-table parity: 2,931 eligible TRB chains and zero mismatches | Empirical interoperability | Indispensable |
| 1E | Ambiguity-aware multi-TRB cell-summary policies | Software behavior | Optional; move to Extended Data if space is limited |

**Panels that should move to Extended Data.** Detailed summary-policy decision
trees, serialization mechanics, subsetting, and concatenation belong in
Extended Data 3. Figure 1 should show only the default ambiguity-aware behavior
in one callout.

**Reviewer objection addressed.** “Is this merely a table wrapper, or does it
behave like a native scverse method without discarding receptor chains?”

**Evidence support.** Yes. The architecture is implemented and tested, and the
public Scirpy example provides exact native/table parity plus H5MU reload
evidence. Panels 1A–1C are explanatory schematics, not performance evidence.

## Figure 2 — Canonical fidelity and scalable deterministic execution

**Single message.** scRFU reproduces the official RFU assignment and scales
deterministically from bounded native objects to a complete 610k-cell
receptor-only atlas, with effective parallel and restart behavior.

| Panel | Content | Evidence role | Priority |
|---|---|---|---|
| 2A | Exact official-RFU parity across ID, label, score, threshold and row reconstruction | Scientific fidelity | Indispensable |
| 2B | Full Wells funnel: 610,429 source cells → 303,088 productive TRB rows → 192,675 unique CDR3 queries; runtime and peak RSS | Full-scale execution | Indispensable |
| 2C | Bounded 1k/10k/25k serial/parallel timing with assignment-hash invariance across chunking and input order | Determinism and parallelism | Indispensable |
| 2D | Fresh versus resumed runtime for Wells 25k and full backend execution | Recovery and cache value | Indispensable |
| 2E | Native Wells 25k/100k/250k runtime, RSS, cached execution and zero table/reload mismatch | Native scale bridge | Indispensable |
| 2F | Baseline, annotated and cell-summary serialized sizes | Storage cost | Optional; preferred in Extended Data 4 |

**Panels that should move to Extended Data.** Detailed chunk-size contrasts,
child-versus-parent RSS, per-stage provenance, and serialized storage overhead
should move to Extended Data 1 and 4. Main panel 2E should retain one memory
trace or point estimate so “native scale” is not presented as runtime alone.

**Reviewer objection addressed.** “Does the portable implementation change
canonical assignments, and can it operate reproducibly at realistic
single-cell scale?”

**Evidence support.** Yes. Exact official parity, full Wells execution, native
25k/100k/250k parity, cache reuse, parallel equivalence, chunk invariance and
order invariance are all frozen. Performance should not be generalized beyond
the audited host.

## Figure 3 — RFU representation properties and frozen feature reuse

**Single message.** RFUs provide a lower-dimensional, less sparse feature
vocabulary than exact receptor identities and that same frozen vocabulary is
reused across independent datasets.

| Panel | Content | Evidence role | Priority |
|---|---|---|---|
| 3A | Exact CDR3, exact clonotype where available, nearest RFU and threshold-qualified RFU feature counts across four public datasets | Representation compression | Indispensable |
| 3B | Distribution of distinct CDR3 sequences per RFU, including medians and singleton fractions | Many-to-one representation structure | Indispensable |
| 3C | Feature count, effective dimensionality and sample-feature sparsity for RFU and comparator representations | Matrix geometry | Indispensable |
| 3D | Frozen-threshold coverage and RFU richness across Wells, GSE190905, GSE157007 and `wu2020_3k` | Reference coverage | Indispensable |
| 3E | Pairwise exact-CDR3 versus nearest/threshold-qualified RFU sharing across datasets | Reusable frozen feature space | Indispensable |
| 3F | Within-dataset RFU sample prevalence distributions | Feature prevalence | Optional; Extended Data |

**Panels that should move to Extended Data.** Full RFU group-size ECDFs,
score quantiles and sample-prevalence distributions may move to Extended Data.
The main compression panel should report absolute counts, not only ratios, to
avoid exaggerating small-cohort effects.

**Reviewer objection addressed.** “Does RFU assignment provide a genuinely
reusable representation, or only relabel exact clonotypes within each cohort?”

**Evidence support.** Yes for the representation-property claim. It does not
show that compression is inherently better, that RFU-sharing implies receptor
equivalence, or that coverage is calibrated OOD detection.

## Figure 4 — Independent transfer, repeated-donor structure and robustness

**Single message.** With the reference and policies frozen, RFU sample vectors
transfer to independent and preregistered held-out cohorts and preserve useful
repeated-donor structure with behavior complementary to Scirpy clonotypes and
conventional representations.

| Panel | Content | Evidence role | Priority |
|---|---|---|---|
| 4A | Wells development/stress, GSE190905 independent repeated-measures, and GSE157007 preregistered held-out roles | Validation design | Indispensable |
| 4B | GSE190905 leave-one-timepoint-out top-1, top-3 and MRR by representation | Donor retrieval | Indispensable |
| 4C | GSE190905 within- versus between-donor cosine by representation | Repeated-donor structure | Indispensable |
| 4D | GSE190905 50%/75% fixed-seed subsampling cosine by representation | Comparative robustness | Indispensable |
| 4E | GSE157007 held-out coverage and nearest/threshold-qualified richness | Frozen transfer | Indispensable |
| 4F | GSE157007 nearest and threshold-qualified sample-vector stability | Held-out robustness | Indispensable |

**Panels that should move to Extended Data.** Jaccard results, all individual
query ranks, diversity components, full comparator matrices, and seed-level
subsampling observations should move to Extended Data 5. Main panels must
retain V/J and Scirpy clonotypes so the complementary result is visible.

**Reviewer objection addressed.** “Does the RFU feature space generalize beyond
the development atlas, and were comparators evaluated on identical samples and
restrictions?”

**Evidence support.** Yes for frozen transfer, technical stability and bounded
two-visit donor representation. RFU does not win every metric: exact CDR3 has
slightly higher MRR, V/J has the highest MRR among the compared count
representations, and coarse V/J/length summaries have higher subsampling
stability partly because they discard receptor detail.

## Figure 5 — Single-cell interpretability and external antigen evidence

**Single message.** Chain-aligned RFUs remain connected to cellular phenotype
and sample context and show group-level coherence of external antigen
annotations beyond size-preserving expectations.

| Panel | Content | Evidence role | Priority |
|---|---|---|---|
| 5A | AIRR chain → RFU annotation → explicit cell/phenotype/sample join → pseudobulk and coupling | Single-cell context architecture | Indispensable schematic |
| 5B | Full Wells RFU-by-cell-type phenotype-coupling heatmap for the 30 most abundant RFUs under nearest and threshold-qualified policies | Descriptive cellular interpretation | Indispensable |
| 5C | Wells phenotype-coupling stability under cell subsampling | Robustness of descriptive linkage | Indispensable |
| 5D | VDJdb exact CDR3 and strict CDR3+V matched-sequence/RFU coverage across three datasets | External evidence coverage | Indispensable |
| 5E | Wells nearest/fractional CDR3 same-antigen pair fraction versus 1,000-permutation size-preserving null | External annotation coherence | Indispensable |
| 5F | Public `wu2020_3k` AIRR-chain → RFU identity → CDR3/CDR3+V query → cell/sample summary, with all 7,544 chains retained | Native evidence linkage | Indispensable |

**Panels that should move to Extended Data.** All 24 VDJdb sensitivity cells,
ambiguity policies, evidence-score distributions, grouping baselines and four
null strata belong in Extended Data 6. The main phenotype heatmap ranks RFUs by
total abundance only, takes the top 30 with ties resolved by RFU label, and
must not select RFUs for visually strong phenotype enrichment.

**Reviewer objection addressed.** “Can RFU annotations be interpreted in the
single-cell ecosystem, and is there any external evidence that sequences
grouped together share annotation structure?”

**Evidence support.** Yes, descriptively. The Wells phenotype analyses retain
cell context without cell-level inferential testing. The VDJdb result supports
external annotation coherence only; database bias, ambiguity and sparse strict
matching remain explicit limitations.

## Extended Data and Supplementary architecture

### Extended Data 1 — Workflow and provenance

- Evidence manifests and artifact/reference hashes across every main evidence
  family.
- Exact-CDR3 deduplication, chunk orchestration and deterministic row
  reconstruction.
- Full run identity, cache compatibility checks and failure/recovery flow.

### Extended Data 2 — Threshold and assignment-policy sensitivity

- Coverage and richness across thresholds, retaining 0.6 as the frozen primary
  threshold.
- Nearest versus threshold-qualified RFU metrics, pseudobulk and phenotype
  coupling.
- Explicit statement that threshold failure is not a calibrated OOD
  probability.

### Extended Data 3 — Scverse object lifecycle

- H5AD/H5MU serialization and reload.
- AnnData/MuData subsetting.
- Compatible concatenation and rejection of incompatible reference/schema
  identities.
- Multi-chain and ambiguity-aware cell-summary behavior.

This is supported by exact tests and native manifests, but there is no separate
frozen plot-ready source table. It should therefore be a technical pass/fail
matrix, not a quantitative performance panel.

### Extended Data 4 — Native memory and storage overhead

- Parent and RFU-child peak RSS at 25k/100k/250k.
- Baseline AIRR, chain-annotated, and optional-summary serialized sizes.
- Serialization and reload times.
- Clear distinction between bounded in-memory native objects and the targeted
  HDF5 route used for full Wells.

### Extended Data 5 — Full comparator results

- Representation dimensionality, sparsity and effective dimension.
- All GSE190905 retrieval metrics and query-rank distributions.
- Cosine and Jaccard within/between summaries.
- Seed-level 50% and 75% downsampling results.
- Explicit display of endpoints where exact CDR3, Scirpy clonotypes, V/J,
  CDR3 length or diversity equal or exceed RFU.

### Extended Data 6 — VDJdb sensitivity and null analyses

- CDR3 versus CDR3+V matching.
- Nearest versus threshold-qualified assignment.
- Fractional versus exclude-ambiguous policies.
- Purity, normalized entropy, same-antigen pair fraction and evidence-score
  distributions.
- RFU, TRBV, CDR3 length, TRBV+length, size-matched random and edit-distance
  grouping comparisons.
- Unrestricted, CDR3-length, TRBV and TRBV+length stratified 1,000-permutation
  nulls, retaining sparse/undefined cells.

### Extended Data 7 — Dependency, CI and interoperability evidence

- Python 3.10/3.11/3.12 Linux matrix.
- Ubuntu, macOS and Windows Python-only smoke.
- Optional dependency, documentation, wheel, sdist and tutorial jobs.
- Current AnnData, MuData, Scirpy, Awkward, NumPy and pandas compatibility
  boundaries.
- Public `wu2020_3k` official-RFU interoperability details.

This is technical-record evidence rather than a scientific result panel.
Official RFU/R execution must remain labelled Linux-validated only.

## Evidence gaps and deliberately excluded panels

No proposed quantitative main panel lacks frozen evidence. The following items
are deliberately not promoted to quantitative panels:

1. **A new UMAP colored by RFU.** The frozen evidence demonstrates native
   context and phenotype linkage, but no separately frozen UMAP-to-RFU source
   table exists. A schematic is sufficient for Figure 5A; a new quantitative
   UMAP should not be fabricated during figure assembly.
2. **Serialization/subsetting/concatenation effect sizes.** These are exact
   software invariants, not biological quantities. Tests and manifests support
   a pass/fail Extended Data matrix.
3. **Deep longitudinal dynamics.** GSE190905 has two visits. GSE345124 QC is
   frozen but RFU execution was explicitly deferred; persistence, expansion,
   contraction and multi-visit trajectory panels are unsupported.
4. **RFU-level antigen specificity.** VDJdb supports annotation coherence and
   null comparisons, not specificity labels for individual RFUs.
5. **Calibrated OOD detection.** Frozen threshold coverage is available, but no
   probability calibration experiment exists.
6. **Universal superiority over Scirpy clonotypes or conventional summaries.**
   The frozen comparator results are complementary and must remain so.

## Is a new experiment required?

**No new experiment is required before manuscript figure design for the stated
central claim.** The five figures can be assembled entirely from the frozen
schematics, source tables and sealed manifests listed in the panel map.

A full GSE345124 analysis would become required only if the scope is expanded
to claim deep longitudinal RFU persistence or multi-visit repertoire dynamics.
It is not required for scalable scverse-native assignment, transferable feature
space, bounded repeated-donor representation, phenotype integration or VDJdb
annotation-coherence claims.
