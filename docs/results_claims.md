# Frozen Results claims

Frozen on 2026-08-24. A claim can move from blocked to demonstrated only after
the corresponding real source table is audited in
[`manuscript_claim_audit.md`](manuscript_claim_audit.md). Software availability
alone is not evidence.

## 1. Deterministic and scalable functional-unit assignment — demonstrated

- **One-sentence claim:** scRFU preserves official standard-mode RFU semantics
  while adding exact sequence deduplication, row reconstruction, bounded
  parallel execution, and resumable chunk processing.
- **Supporting dataset:** adversarial official-parity fixture and Wells 1k/10k/25k.
- **Figure panel:** Figure 1B–1F.
- **Source table:** `original_rfu_parity/figure1_original_rfu_parity.tsv` and
  `wells_bounded/figure1_scaling_source.tsv`.
- **Quantitative result:** 10/10 eligible parity rows matched with zero score,
  threshold, and reconstruction discrepancy; 25k 5,000-query chunks took
  45.650 s with one worker versus 25.859 s with two, while resume took
  0.209–0.222 s and all 25k configurations shared one assignment hash.
- **Uncertainty/test:** exact equivalence at `1e-12`; runtime is a single
  hardware-specific observation, not an inferential estimate.
- **Comparator:** official upstream `AssignRFUs()` and serial/fresh execution.
- **Caveat:** the parity fixture is bounded and scaling is not extrapolated to
  the full atlas.

## 2. Frozen RFU assignment transfers across independent cohorts — demonstrated for technical transfer

- **One-sentence claim:** the unchanged official RFU reference produced usable
  and similar threshold coverage across Wells, independent GSE190905, and
  preregistered held-out GSE157007.
- **Supporting dataset:** Wells 25k, GSE190905, GSE157007.
- **Figure panel:** Figure 3A–3C.
- **Source table:** Wells assignment-policy summary, GSE190905 run manifest, and
  GSE157007 transfer summary.
- **Quantitative result:** threshold coverage was 76.94%, 78.31%, and 76.88%;
  nearest RFU richness was 3,963, 4,465, and 4,898 respectively.
- **Uncertainty/test:** descriptive cohort-specific estimates; no cross-cohort
  equality test.
- **Comparator:** nearest versus threshold-qualified assignment under the same
  frozen reference.
- **Caveat:** similar coverage is not proof of biological equivalence,
  population generalization, or out-of-distribution validity.

## 3. RFU representations remain stable under substantial repertoire subsampling — demonstrated

- **One-sentence claim:** RFU abundance and held-out sample vectors degrade
  gradually rather than discontinuously under deterministic bounded
  subsampling.
- **Supporting dataset:** Wells 25k and preregistered GSE157007.
- **Figure panel:** Figure 1H and Figure 3F.
- **Source table:** `wells_bounded/figure1_robustness_source.tsv` and
  `public_data/GSE157007/heldout_validation/subsampling_stability.tsv`.
- **Quantitative result:** Wells nearest-RFU median cosine was
  0.526/0.730/0.885 at 25/50/75% cell retention; held-out nearest-RFU mean
  sample cosine was 0.936 at 50% and 0.976 at 75%.
- **Uncertainty/test:** individual values from three fixed seeds with medians or
  means as declared by each source table.
- **Comparator:** full bounded result and conventional held-out representations.
- **Caveat:** multinomial abundance resampling is not called physical
  sequencing-depth simulation; bounded subsampling does not establish clinical
  robustness.

## 4. scRFU links receptor functional units with single-cell phenotypes — demonstrated descriptively

- **One-sentence claim:** the bounded Wells workflow connects RFU assignments
  with cell-type annotations and retains much of that coupling under
  deterministic subsampling.
- **Supporting dataset:** Wells 25k.
- **Figure panel:** Figure 3G–3H.
- **Source table:** nearest and threshold-qualified Wells phenotype-coupling
  tables and `wells_bounded/phenotype_coupling_stability.tsv`.
- **Quantitative result:** the receptor-bearing subset covered 16 cell types,
  21 donors, and 209 libraries; nearest-policy coupling cosine was
  0.506/0.708/0.863 at 25/50/75% cell retention.
- **Uncertainty/test:** three fixed seeds; descriptive coupling profiles without
  cell-level p-values.
- **Comparator:** nearest versus threshold-qualified RFU policy plus exact CDR3,
  V, length, and diversity summaries where available.
- **Caveat:** the bounded subset is not claimed biologically representative and
  coupling does not imply antigen specificity or phenotype causality.

## 5. RFU representations resolve governed longitudinal donor structure and repertoire dynamics — blocked

- **One-sentence claim:** `[BLOCKED: LONGITUDINAL RESULT]`.
- **Supporting dataset:** governed six-volunteer repeated-measures cohort.
- **Figure panel:** Figure 2, currently blocked in full.
- **Source table:** none.
- **Quantitative result:** none.
- **Uncertainty/test:** prespecified donor-aware exact permutation, donor
  bootstrap, and leave-one-donor-out sensitivity have not run.
- **Comparator:** RFU, exact CDR3, supplied clonotype, V, J, length, diversity,
  and feasible edit-distance groups with identical candidate sets.
- **Caveat:** explicit runtime paths are unset; six participants could support a
  deep methods demonstration only, never population-level aging inference.

The paired GSE190905 analysis may be reported separately as independent
technical evidence: nearest-RFU mean cosine was 0.491 within patient versus
0.063 between patients, and cosine retrieval was top-1 0.667, top-3 0.833, MRR
0.764 across 12 queries. It does not substitute for the governed cohort.

## 6. RFUs show antigen-evidence coherence beyond simple grouping baselines — blocked

- **One-sentence claim:** `[BLOCKED: VDJDB RESULT]`.
- **Supporting dataset:** real release-pinned VDJdb matched to Wells,
  GSE190905, and GSE157007.
- **Figure panel:** Figure 4 primary design, currently blocked.
- **Source table:** none.
- **Quantitative result:** none.
- **Uncertainty/test:** at least 1,000 size-preserving unrestricted,
  CDR3-length-stratified, TRBV-stratified, and feasible joint-stratified
  permutations remain pending.
- **Comparator:** TRBV, CDR3 length, TRBV+length, size-matched random partitions,
  and feasible edit-distance clusters.
- **Caveat:** `VDJDB_PATH` and `VDJDB_RELEASE` are unset. Even a positive result
  would support antigen-label coherence among distinct sequences, not claim
  that an RFU is antigen-specific.
