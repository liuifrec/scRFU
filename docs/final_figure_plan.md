# Final manuscript figure plan

Frozen on 2026-08-24 from evidence that has actually been executed. Runtime
paths are relative to `SCRFU_MONTH2_OUTDIR`. Panel lettering may change during
layout, but panel definitions and claims may not be expanded without updating
the claim audit and source-data inventory.

## Figure 1 — Framework and technical validation

| Panel | Source table | x / y or displayed variable | Unit | Uncertainty | Comparator | Caption claim | Required caveat |
|---|---|---|---|---|---|---|---|
| 1A workflow | `docs/methods_freeze.md`; no quantitative source | canonical receptor input → exact-CDR3 deduplication → official frozen RFU assignment → row reconstruction → sample/phenotype summaries | method component | none | none | scRFU applies an immutable upstream RFU reference while preserving receptor-row identity. | Schematic is not evidence of performance. |
| 1B official parity | `original_rfu_parity/figure1_original_rfu_parity.tsv` | comparison category / exact matches or mismatch count; inset maximum score error | receptor row | exact result; tolerance `1e-12` | official upstream `AssignRFUs()` | Canonical scRFU standard mode reproduced official assignments for every eligible adversarial fixture row. | Bounded fixture, not repertoire-scale accuracy. |
| 1C deduplication | `wells_bounded/figure1_scaling_source.tsv` | source-cell subset / receptor rows and unique queries, or unique-query/input ratio | receptor row and unique CDR3 | deterministic point estimates | raw eligible rows | Exact-CDR3 deduplication reduced backend queries while retaining row-level reconstruction. | Benefit depends on repertoire duplication. |
| 1D serial/parallel scaling | same scaling table | source cells or unique queries / wall time; color worker count; facet chunk size | one execution | observed run only; no inferential interval | one versus two process workers | Two-worker chunk execution reduced bounded 25k wall time with identical assignments. | Hardware-specific; no full-atlas extrapolation. |
| 1E restart and invariance | same scaling table | configuration / wall time or equality indicator | one execution/configuration | exact assignment SHA256 | fresh versus resumed; serial versus parallel; chunk/order variants | Resume reused completed chunks and every meaningful execution configuration produced the same 25k assignment hash. | Resume requires unchanged run identity and artifacts. |
| 1F memory | same scaling table | source cells/configuration / Python and backend peak RSS | execution process | observed peaks | serial/parallel configurations | Bounded memory use was measured for Python and the RFU backend. | Peaks are platform-specific and are not additive. |
| 1G reference coverage | `wells_bounded/threshold_sensitivity.tsv`; `wells_bounded/single_cell_analysis/assignment_policy_summary.tsv` | threshold / pass fraction; policy / assigned cells and RFUs | receptor-bearing cell | descriptive point estimates | nearest versus threshold-qualified assignment | Nearest and threshold-qualified semantics are reported separately, with 76.94% Wells 25k coverage at the frozen threshold. | Threshold sensitivity is descriptive; 0.6 remains frozen. |
| 1H robustness | `wells_bounded/figure1_robustness_source.tsv` | retained fraction / cosine, Spearman, Jaccard, top-k overlap, MAE, or rank displacement | RFU abundance vector | three fixed seeds; show individual seeds and median | full Wells 25k bounded reference | RFU abundance patterns degraded gradually under deterministic cell and sequence subsampling. | Multinomial abundance resampling is not physical sequencing-depth simulation. |

## Figure 2 — Longitudinal validation

**Status: blocked.** Figure 2 is not frozen as an evidence figure because the
three required governed runtime paths are unset. No synthetic result may occupy
these panels.

If the governed run becomes available, the prespecified panels are cohort QC,
within/between donor similarity, leave-one-timepoint-out donor retrieval,
RFU persistence/dynamics, compartment trajectories, technical resampling, and
donor-level leave-one-out/bootstrap/exact permutation. Each panel must use a
biological sample, donor, or donor–RFU trajectory as its unit; never a cell as an
inferential replicate.

## Figure 3 — Cross-cohort transfer and single-cell interpretation

| Panel | Source table | x / y or displayed variable | Unit | Uncertainty | Comparator | Caption claim | Required caveat |
|---|---|---|---|---|---|---|---|
| 3A cohort roles/design | `docs/cross_dataset_evidence_summary.md`; public acquisition manifests | dataset / design, receptor rows, unique CDR3s | cohort | none | none | Wells, GSE190905, and GSE157007 provide distinct development, independent-validation, and held-out roles. | Roles cannot be blurred; cohort purposes and platforms differ. |
| 3B frozen-reference coverage | Wells assignment-policy summary; `public_data/GSE190905/validation/run_manifest.json`; `public_data/GSE157007/heldout_validation/transfer_summary.tsv` | dataset / threshold-pass fraction | receptor-bearing cell within dataset | descriptive point estimate | nearest versus threshold-qualified semantics | The unchanged official reference yielded similar bounded threshold coverage across three public cohorts. | Similar coverage does not establish biological equivalence or out-of-distribution validity. |
| 3C RFU richness | same sources plus cohort RFU summaries | dataset / nearest and threshold-qualified RFUs | cohort | descriptive point estimate | assignment policy | Thousands of RFUs were observed in each transferred cohort under the same reference. | Richness depends on cohort size and sampling depth; do not compare as an outcome without normalization. |
| 3D independent paired retrieval | `cross_cohort/gse190905_donor_retrieval_comparators.tsv` | representation / top-1, top-3, MRR, correct-donor rank | patient-time query | 12 leave-one-timepoint-out queries; show all query ranks where space permits | exact CDR3, clonotype, V, J, length, diversity | RFU representations retained donor-retrieval information with performance complementary to conventional representations. | Six treated patients; descriptive methods validation, not population inference. |
| 3E independent within/between structure | `cross_cohort/gse190905_within_between_comparators.tsv` | donor relation / cosine or distance | patient-time sample pair | show six within-donor and 60 between-donor pairs; no naive pairwise p-value | same representations | Within-patient RFU cosine similarity exceeded the between-patient mean in this paired cohort. | Pairs are dependent; no donor-aware inferential test was run. |
| 3F held-out stability | `public_data/GSE157007/heldout_validation/subsampling_stability.tsv` and summary | retained fraction / sample-vector cosine | biological sample | three seeds across 17 samples; show distribution and mean | RFU nearest/threshold, exact CDR3, clonotype, V, J, length, diversity | Preregistered held-out RFU vectors remained stable at 50% and 75% deterministic subsampling. | Cross-sectional one-sample-per-donor design; no retrieval endpoint. |
| 3G Wells phenotype coupling | `wells_bounded/single_cell_analysis/nearest/rfu_phenotype_coupling.tsv`; threshold-qualified counterpart | cell type / RFU coupling, specificity, or abundance | RFU×cell-type profile | descriptive; donor/library summaries where used | nearest versus threshold-qualified; conventional repertoire summaries | scRFU connects receptor functional-unit assignments with single-cell phenotype annotations at bounded atlas scale. | No cell-level inferential p-values and no claim of subset representativeness. |
| 3H phenotype-coupling robustness | `wells_bounded/phenotype_coupling_stability.tsv` | retained fraction / coupling cosine and dominant-phenotype agreement | RFU×phenotype profile | three fixed seeds | nearest versus threshold-qualified | Phenotype-coupling profiles were progressively stable under cell and sequence subsampling. | Descriptive stability, not proof of phenotype causality. |

## Figure 4 — Antigen-evidence coherence or fallback

**Primary status: blocked.** A real release-pinned VDJdb reference is not
configured. Figure 4 may become an antigen-evidence coherence figure only if
the audited match coverage supports interpretable group-level comparisons. It
must not describe any RFU as antigen-specific.

The prespecified primary panels are matched-sequence coverage, RFUs with at
least two independently matched sequences, purity/normalized entropy/
same-antigen-pair fraction, ambiguity-policy sensitivity, size-preserving and
receptor-property-stratified permutation nulls, and comparisons with TRBV,
length, TRBV+length, size-matched random groups, and feasible edit-distance
groups. The unit is a distinct `unique_sequence_antigen` evidence record or RFU;
uncertainty is the empirical distribution from at least 1,000 permutations.

If real VDJdb coverage is insufficient, Figure 4 becomes an expanded technical
robustness/transfer figure using only the already audited Wells robustness,
GSE190905 paired retrieval, GSE157007 held-out stability, and assignment-policy
sensitivity tables. BCR is not a fallback and is deliberately excluded from the
first manuscript.
