# Archived 2026-08-24 manuscript: frozen-reference receptor states

Historical draft retained unchanged below this notice. The active manuscript is
[`../scRFU_methods_manuscript.md`](../scRFU_methods_manuscript.md). Its cohort
boundary and results supersede this draft; old outcome-dependent selections are
not authorized for the radiation-methods workstream.

Draft status: evidence-frozen manuscript skeleton, 2026-08-24. Target format:
Cell Reports Methods; prose remains journal-neutral for later adaptation.

## Highlights

- scRFU exactly reproduces canonical standard-mode assignments from the public RFU method.
- Exact-sequence deduplication, bounded parallelism, and resumable execution preserve assignments.
- One frozen T-cell RFU reference transfers across three distinct public receptor datasets.
- RFU assignments connect receptor states with single-cell phenotype annotations.

## Summary

Comparing adaptive immune repertoires across samples is difficult because exact
receptor sequences are sparse, whereas broad gene-usage and diversity summaries
discard sequence-level structure. Receptor functional units (RFUs) provide an
intermediate representation by assigning receptor sequences to a fixed
reference of sequence-defined states. We developed scRFU, a reproducible Python
workflow that normalizes single-cell receptor inputs, deduplicates identical
queries, delegates canonical scoring to the public RFU implementation, and
reconstructs assignments to the original cell order. In an adversarial parity
test, all 10 eligible rows reproduced the official RFU identifier, label,
threshold status, order, and score, with maximum score difference zero. In a
bounded 25,000-cell Wells atlas subset, 12,415 productive TRB records reduced to
10,038 unique CDR3 queries; two-worker execution reduced wall time from 45.65 to
25.86 s for 5,000-query chunks, and cached resumption required approximately
0.21 s while preserving an identical assignment hash. The unchanged RFU
reference yielded threshold coverage of 76.94%, 78.31%, and 76.88% in Wells,
GSE190905, and preregistered held-out GSE157007, respectively. RFU abundance and
single-cell phenotype-coupling profiles degraded gradually under deterministic
subsampling. Conventional representations showed complementary rather than
uniformly inferior behavior. These results establish scRFU as a deterministic,
bounded, transferable T-cell receptor-state workflow. Governed longitudinal and
release-pinned antigen-evidence analyses remain explicitly blocked pending their
runtime inputs.

## Introduction

Adaptive immune repertoires contain many distinct receptor sequences, most of
which are rare or private to an individual. Exact clonotype comparisons retain
fine sequence information but often produce sparse sample-by-feature matrices.
Conversely, V/J-gene usage, CDR3 length, and scalar diversity summaries are
portable but compress potentially meaningful sequence variation. A useful
methods representation should bridge these levels while remaining stable,
interpretable, and applicable without refitting on every validation cohort.

The public RFU method defines a reference-based representation for T-cell
receptor beta-chain CDR3 amino-acid sequences [RFU METHOD REF]. Its scientific
value for single-cell and cross-cohort work depends not only on the assignment
algorithm but also on reproducible data normalization, exact handling of
duplicate sequences, source-order reconstruction, bounded execution, transparent
threshold semantics, and frozen-reference validation. These requirements are
particularly important when millions of cell records contain fewer unique
receptor queries and when a validation cohort must not influence the reference.

We built scRFU to provide this reproducible analysis layer without refitting or
bundling the upstream RFU reference. We evaluated canonical parity, bounded
runtime and memory, chunk/parallel/cache invariance, deterministic robustness,
single-cell phenotype linkage, independent cohort transfer, preregistered
held-out stability, and conventional repertoire comparators. The analyses focus
on technical transferability and descriptive biological structure. They do not
assume that RFUs outperform every representation or that an RFU denotes one
antigen specificity.

## Results

### scRFU provides deterministic and scalable functional-unit assignment

scRFU converts receptor records to canonical productive primary-TRB queries,
deduplicates identical eligible CDR3 amino-acid sequences, executes the
unchanged public RFU standard-mode scorer, and reconstructs assignments to the
source records in stable order (Figure 1A). The frozen reference was the public
RFU commit `ff7a1ea6aca28444a7ff8ec7ef09745e7341f99d`, with its three required
artifacts identified by SHA256 in the Methods.

The parity input contained unique and duplicated CDR3s, identical CDR3s paired
with different TRBV values, non-C-starting sequences, below-threshold queries,
repeated cells, and shuffled order. Of 12 input rows, 10 were eligible and
represented five unique eligible CDR3s. Canonical scRFU and the official
implementation agreed for all 10 eligible RFU identifiers and labels, with zero
threshold, reconstruction, or order mismatch and maximum RFU-score difference
zero at tolerance `1e-12` (Figure 1B).

Deterministic Wells subsets of 1,000, 10,000, and 25,000 source cells contained
499, 4,861, and 12,415 productive TRB records, which reduced to 453, 4,130, and
10,038 unique queries (Figure 1C). At 25,000 cells and 5,000-query chunks, one
worker completed RFU execution in 45.650 s and two workers in 25.859 s. Resumed
runs required 0.209–0.222 s and reused every completed chunk. Serial, parallel,
chunk-size, shuffled-input, and resume configurations reconstructed the same
assignment SHA256 (Figure 1D–F). Runtime and memory values describe this bounded
machine-specific execution and are not full-atlas extrapolations.

### Frozen RFU assignment transfers across independent receptor cohorts

We used the same reference, standard mode, and 0.6 threshold for the bounded
Wells atlas, independent GSE190905 cohort, and preregistered held-out GSE157007
cohort. These datasets contained 12,415, 27,655, and 60,125 receptor rows and
10,038, 16,391, and 43,082 unique CDR3s, respectively. Threshold-qualified
coverage was 76.94%, 78.31%, and 76.88%, while nearest assignment observed
3,963, 4,465, and 4,898 RFUs (Figure 3A–C). Similar coverage under an unchanged
reference supports technical transfer; it does not demonstrate biological
equivalence among cohorts.

GSE190905 provided 12 paired samples from six treated patients. Across 12
leave-one-timepoint-out queries, nearest-RFU cosine retrieval achieved top-1
accuracy 0.667, top-3 accuracy 0.833, and mean reciprocal rank 0.764. Mean
nearest-RFU cosine similarity was 0.491 for the six within-patient pairs and
0.063 for 60 between-patient pairs (Figure 3D–E). These dependent descriptive
pairs were not subjected to a naive pair-level significance test. Exact CDR3,
clonotype, V, J, length, and diversity representations sometimes matched or
exceeded individual RFU endpoints, supporting complementary rather than
universal-superiority interpretation.

### RFU representations remain stable under substantial repertoire subsampling

In Wells 25k, median nearest-RFU abundance-vector cosine similarity to the full
bounded reference was 0.526, 0.730, and 0.885 at 25%, 50%, and 75% cell
retention; corresponding median Spearman correlations were 0.510, 0.713, and
0.873 (Figure 1H). In preregistered held-out GSE157007, the mean nearest-RFU
sample-vector cosine was 0.936 at 50% and 0.976 at 75% deterministic
subsampling; threshold-qualified values were 0.922 and 0.970 (Figure 3F). Each
summary used three fixed seeds. These perturbations assess bounded stability;
multinomial abundance resampling is not interpreted as physical sequencing-depth
simulation without appropriate read or UMI information.

### scRFU links receptor functional units with single-cell immune phenotypes

The Wells 25k receptor-bearing subset spanned 16 cell types, 21 donors, and 209
libraries. Joining only required observation metadata enabled RFU abundance,
richness, convergence, multiplicity, donor prevalence, pseudobulk, and
RFU-by-cell-type coupling summaries without loading an expression matrix. Under
nearest assignment, the coupling-profile cosine relative to the full bounded
reference was 0.506, 0.708, and 0.863 at 25%, 50%, and 75% cell retention
(Figure 3G–H). Nearest and threshold-qualified analyses are presented
separately. These descriptive associations connect receptor states with cell
annotations but do not establish phenotype causality, antigen specificity, or
population representativeness.

### RFU representations resolve longitudinal donor structure and repertoire dynamics

`[BLOCKED: LONGITUDINAL RESULT]`

The governed six-volunteer input paths are not configured. No cohort count,
within/between-donor test, retrieval statistic, persistence category,
compartment contrast, bootstrap, or permutation result is reported. The paired
GSE190905 analysis is retained as independent technical evidence and is not a
substitute for the prespecified deep repeated-measures cohort.

### RFUs show antigen-evidence coherence beyond simple receptor groupings

`[BLOCKED: VDJDB RESULT]`

No release-pinned VDJdb runtime reference is configured. No match rate,
coherence statistic, evidence-score sensitivity, or permutation null is
reported. If this analysis becomes available, its claim is limited to antigen
label coherence among distinct sequences assigned to the same RFU, not RFU
antigen specificity.

## Discussion

The completed analyses support three aspects of the scRFU methods claim. First,
the workflow retains canonical upstream behavior while adding deterministic
deduplication, reconstruction, parallel execution, and cache reuse. Second, a
single frozen reference yielded substantial threshold coverage across three
public datasets with distinct roles, including a preregistered held-out cohort.
Third, the representation can be joined to single-cell phenotype metadata and
remains progressively stable under bounded subsampling.

RFU and conventional repertoire representations summarize different
information. Exact CDR3s preserve identity but are sparse, gene usage and CDR3
length are compact, and RFUs aggregate sequences through a fixed reference.
Observed comparator results therefore motivate a complementary-methods claim,
not a ranking in which RFU must dominate every task. A frozen representation is
most useful when portability, reproducibility, and intermediate granularity are
the target.

The central evidence presently supports technical transfer and descriptive
single-cell interpretation. A governed repeated-measures run is still required
before making the manuscript's proposed longitudinal claim. Release-pinned
VDJdb analysis is optional for the central transfer story but required for any
antigen-coherence result. If its coverage is insufficient, Figure 4 will use
prespecified robustness and transfer evidence instead.

## Methods

### Study design and dataset roles

Wells was used as a bounded development and single-cell demonstration;
GSE190905 as an independent paired technical validation; and GSE157007 as a
preregistered held-out cross-sectional validation. Dataset roles were fixed
before synthesis. The governed longitudinal cohort remains unexamined because
its explicit runtime inputs are absent.

### RFU reference and artifact integrity

All analyses used public RFU commit
`ff7a1ea6aca28444a7ff8ec7ef09745e7341f99d`. SHA256 values were
`92c2faa33f2e7f60d6470ad7dfd653eb0ca54859aa198ed53ad0c876aae8640b`
for `RFU.R`,
`820fb71428913974e994543cbbbaa591f54355cdb73d651f068f0e691a6cdffd`
for the trimer reference, and
`64783074360edeca84b3f49ab9d682263e954752fc09bb208edaabde0fd82553`
for the centroid/reference. Assets were used from an external untouched
checkout and were not bundled with scRFU.

### Receptor normalization, eligibility, and reconstruction

Analyses selected productive primary TRB records under each frozen dataset
manifest. CDR3 amino-acid strings were whitespace-normalized but not
approximately collapsed. Exact eligible CDR3s were queried once and results
were reconstructed to all source records in stable order. TRBV was not part of
the standard-mode query key. Official standard-mode eligibility and the 0.6
threshold were retained. Nearest assignment reports the nearest reference state
regardless of threshold; threshold-qualified assignment retains only scores
passing 0.6.

### Chunking, parallelism, cache reuse, and provenance

Chunking occurred after deduplication. Process workers executed independent
chunks, and reconstruction followed query order rather than task-completion
order. Cache identity included ordered query identifiers and sequences,
artifact and wrapper hashes, RFU mode, threshold, and execution parameters.
Every source table recorded input identifiers or hashes, the frozen reference,
software version, parameters, seed where applicable, and date.

### Bounded Wells sampling and single-cell summaries

Cells were selected by deterministic hash-ranked random sampling with seed
20260824 at 1,000, 10,000, and 25,000 cells. A proportional stratified 25,000
subset was evaluated as a sampling diagnostic but was not substituted after
viewing results. Observation and receptor metadata were selectively read;
expression matrices were not loaded. Single-cell summaries were calculated on
reconstructed receptor-bearing cells, without cell-level inferential p-values.

### Robustness analyses

Cell and sequence subsampling and multinomial abundance resampling used
fractions 0.25, 0.50, 0.75, and 1.00 and seeds 20260824–20260826. Perturbed
vectors were compared with the full bounded reference using prespecified
Spearman, cosine, Jaccard, top-k overlap, absolute-error, rank-displacement, and
coverage metrics. Threshold, chunk, processing-order, parallel, and resume
sensitivity were evaluated without changing the frozen primary threshold.

### Conventional repertoire comparators

Comparator matrices used identical samples and receptor rows for exact CDR3,
clonotype where supplied, normalized V usage, normalized J usage, and CDR3
length. Shannon and Simpson diversity were computed from the same abundance
inputs. Edit-distance connected components at amino-acid distance at most one
were prespecified only for at most 2,000 unique sequences and were recorded as
unavailable above this frozen computational bound.

### Independent and held-out validation

GSE190905 was evaluated with the unchanged RFU reference and paired sample
metadata. Leave-one-timepoint-out retrieval used identical candidate sets for
all representations. GSE157007 RFU reference, threshold, chain, weighting,
subsampling, metrics, comparators, and metadata harmonization were registered in
`docs/heldout_gse157007_preregistration.json` before held-out evaluation and
were not tuned afterward.

### Longitudinal analyses

`[BLOCKED: LONGITUDINAL RESULT]` Primary proportion normalization, sample-pair
metrics, leave-one-timepoint-out retrieval, count-based dynamics categories,
compartment comparisons, resampling, donor bootstrap, leave-one-donor-out, and
donor-label permutation are prespecified in `docs/methods_freeze.md`. The donor,
not the cell or sequence, is the biological replicate.

### Antigen-evidence matching and null models

`[BLOCKED: VDJDB RESULT]` Exact CDR3 and exact CDR3-plus-normalized-V matching,
`unique_sequence_antigen` aggregation, fractional and exclude-ambiguous
policies, evidence-score sensitivity, and four size-preserving 1,000-permutation
nulls are prespecified in `docs/methods_freeze.md`.

### Software and statistical reporting

Completed source tables were generated with scRFU 0.1.0 on 2026-08-24.
Technical equivalence uses exact assignment comparison and score tolerance
`1e-12`. Descriptive biological summaries do not receive cell- or sequence-level
p-values. Hardware-specific runtime and RSS are reported as observed values.

## Limitations

The parity fixture is adversarial but bounded, and the atlas benchmark stops at
25,000 source cells. Similar threshold coverage across cohorts does not prove
biological equivalence or universal out-of-distribution performance. The Wells
subset is not claimed representative of the full atlas. GSE190905 has six
patients and supports a paired technical demonstration rather than population
inference; GSE157007 is cross-sectional with one sample per donor. Comparator
performance is task dependent, and edit-distance grouping is unavailable above
the frozen quadratic-work bound. Real governed longitudinal and VDJdb analyses
remain absent. No antigen-specificity claim is supported. This first manuscript
is limited to T-cell receptor beta-chain RFU analysis.

## Data availability

GSE190905 and GSE157007 are publicly available through NCBI GEO under their
stated accessions. Public dataset URLs, expected files, hashes, and acquisition
instructions are recorded in repository manifests. The Wells atlas is a public
source dataset used through a user-supplied local copy. Public source data are
not redistributed by scRFU. Governed longitudinal inputs, if analyzed, will
remain controlled; only disclosure-reviewed aggregate derived tables without
participant identifiers may be shared. VDJdb is not bundled or redistributed.

## Code availability

Source code is developed publicly at `https://github.com/liuifrec/scRFU` under
the repository license. The current evidence run uses version 0.1.0; a coherent
0.4.0 release, archive, and DOI are planned but do not yet exist. RFU assets are
obtained separately from the public upstream repository.

## Figure legends

### Figure 1. scRFU workflow and technical validation

(A) Frozen-reference workflow. (B) Exact standard-mode parity with official RFU
for the adversarial fixture. (C) Receptor-row and unique-query counts in bounded
Wells subsets. (D) Serial and two-worker timing. (E) resume and execution-order
invariance. (F) Python and backend peak RSS. (G) nearest and
threshold-qualified coverage. (H) deterministic robustness. Values are bounded
technical measurements; runtime and memory are platform specific.

### Figure 2. Governed longitudinal methodological validation

`[BLOCKED: LONGITUDINAL RESULT]` No evidence figure is frozen until explicit
governed inputs have been validated and the prespecified donor-level analyses
have run.

### Figure 3. Frozen-reference transfer and single-cell interpretation

(A) Cohort roles and design. (B) threshold-qualified coverage under one
reference. (C) RFU richness by assignment policy. (D) paired GSE190905 donor
retrieval and conventional comparators. (E) within- and between-patient
similarity. (F) preregistered held-out GSE157007 subsampling stability. (G) Wells
RFU–cell-type coupling. (H) coupling robustness. Cohort-specific biological
designs are not pooled as one population.

### Figure 4. Antigen-evidence coherence or prespecified fallback

`[BLOCKED: VDJDB RESULT]` If release-pinned VDJdb coverage is adequate, panels
will report matching, coherence, ambiguity sensitivity, permutation nulls, and
simple grouping comparators without antigen-specificity language. Otherwise the
figure will contain only already audited expanded robustness and transfer
analyses.
