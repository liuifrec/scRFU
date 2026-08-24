# Methods freeze

Freeze date: 2026-08-24. This document fixes the primary analysis definitions
before either governed longitudinal data or a release-pinned VDJdb reference is
examined. Later operational bug fixes may be recorded, but primary definitions
must not change in response to results. All runtime outputs remain outside Git.

## Frozen RFU reference and receptor processing

- **Reference:** untouched public RFU commit
  `ff7a1ea6aca28444a7ff8ec7ef09745e7341f99d` from
  `https://github.com/s175573/RFU`.
- **Artifacts:** `RFU.R` SHA256
  `92c2faa33f2e7f60d6470ad7dfd653eb0ca54859aa198ed53ad0c876aae8640b`;
  trimer SHA256
  `820fb71428913974e994543cbbbaa591f54355cdb73d651f068f0e691a6cdffd`;
  centroid/reference SHA256
  `64783074360edeca84b3f49ab9d682263e954752fc09bb208edaabde0fd82553`.
- **Frozen-reference ID:**
  `scrfu-ref-db43f0ffe17d9fe6ec4d0d36231e08b48e52f5b508b7039cc08024fdf0158518`.
- **Chain:** productive primary TRB. When several productive TRB records map to
  one cell, the dataset-specific deterministic primary-chain rule recorded in
  its acquisition or run manifest is used; the rule is not changed per cohort.
- **Normalization:** trim surrounding whitespace, normalize missing values, and
  retain the source row and source identifier for reconstruction. CDR3 amino-acid
  text is not error-corrected or approximately collapsed. V-gene alleles are
  stripped only for analyses that explicitly use normalized V genes.
- **Eligibility:** the official standard-mode criterion is retained, including
  the canonical C-start requirement. Ineligible rows remain reconstructable but
  do not receive an RFU result.
- **Deduplication:** exact eligible CDR3 amino-acid strings are the query key.
  Identical CDR3s are executed once and reconstructed to every source row in
  stable input order. TRBV is not part of the standard-mode deduplication key
  because it is not an active input to the official standard-mode encoder.
- **Assignment:** standard mode and score threshold 0.6. Nearest assignment
  always reports the nearest RFU and score; threshold-qualified assignment
  retains that assignment only when the score passes 0.6. The threshold does
  not trigger reassignment.
- **Numerical equivalence:** IDs, labels, eligibility, threshold status, and
  reconstruction must be exact; score comparisons use absolute tolerance
  `1e-12`.

## Execution, coverage, and cache semantics

- Query chunking operates only after exact CDR3 deduplication. The primary large
  bounded configuration uses 5,000 unique queries per chunk and two process
  workers; serial and other meaningful chunk sizes are equivalence checks.
- Parallel execution uses a process executor. Output is reconstructed by stable
  query and source order, never completion order.
- A cached chunk is reusable only when its ordered identifiers and sequences,
  RFU artifact hashes, wrapper hash, mode, threshold, and execution parameters
  match. Resume must produce the same assignment hash as fresh execution.
- Reference coverage is the number of eligible receptor rows whose nearest
  score passes 0.6 divided by all eligible receptor rows. Nearest-RFU richness
  and threshold-qualified RFU richness are reported separately.
- Runtime and peak RSS are descriptive for the recorded machine and bounded
  input; they are not extrapolated to the full atlas.

## Single-cell and cross-cohort summaries

- Receptor-bearing cells are the counting unit for cell-weighted RFU abundance.
  Exact CDR3 richness, multiplicity, convergence, donor prevalence, phenotype
  coupling, and pseudobulk are computed after receptor-to-cell reconstruction.
- Nearest and threshold-qualified policies are always reported separately.
- Only interpretation metadata are read from H5AD. Expression matrices and
  `raw/X` are not loaded for these analyses.
- Cross-cohort comparisons are restricted to explicitly harmonized technical
  metrics. Cohort roles remain development demonstration (Wells), independent
  validation (GSE190905), held-out validation (GSE157007), and governed
  longitudinal demonstration (pending).

## Robustness freeze

- Cell and sequence subsampling fractions are 0.25, 0.50, 0.75, and 1.00.
- Primary fixed seeds are 20260824, 20260825, and 20260826.
- Abundance resampling is multinomial resampling of observed abundance unless
  an input contains validated physical read/UMI depth; it is not called physical
  sequencing-depth simulation otherwise.
- Each perturbation is compared with its full bounded reference using Spearman
  correlation, cosine similarity, Jaccard overlap, top-k overlap, mean absolute
  abundance error, rank displacement, and reference-coverage change. Phenotype
  coupling stability is included only where cell metadata support it.
- Chunk-size, processing-order, serial/parallel, and resume invariance require
  exact assignments and score agreement at `1e-12`.

## Governed longitudinal analysis freeze

These choices are prespecified but unexecuted while the three required runtime
paths are unset.

- The biological replicate is the donor. Individual cells, sequences, samples,
  and timepoints are not independent biological replicates.
- Primary sample representations are cell-weighted feature proportions. Raw
  counts are retained for cohort QC and the count-based persistence classifier.
  Nearest assignment is primary and threshold-qualified assignment is a
  required sensitivity analysis.
- Same-donor across-time and different-donor sample pairs are compared with
  cosine similarity, binary Jaccard, weighted Jaccard, Bray-Curtis dissimilarity,
  and Jensen-Shannon distance.
- Donor retrieval is leave-one-timepoint-out with identical candidates for all
  representations. Cosine is primary, Jaccard is secondary, `k=3` is the
  prespecified top-k endpoint, and top-1 accuracy, top-k accuracy, mean
  reciprocal rank, and correct-donor rank are reported.
- Comparator representations are exact CDR3, supplied clonotype when present,
  normalized V-gene usage, normalized J-gene usage, CDR3-length usage,
  conventional diversity summaries, and edit-distance groups only within the
  frozen feasible bound described below.
- Persistence/dynamics are classified on counts with at least two observed
  timepoints, persistence in at least 50% of visits, appearance/disappearance
  threshold zero, twofold expansion/contraction, and pseudocount 0.5. Primary
  minimum abundance is one receptor-bearing cell; sensitivity uses two and five
  cells. Categories are persistent stable, persistent variable, expanding,
  contracting, appearing, disappearing, intermittent, and insufficient
  coverage.
- Compartment analyses use the configured generic compartment field. CD4/CD8
  labels are used in tables only if those values are present. Primary analyses
  include within-compartment stability, cross-compartment similarity, shared
  persistent RFUs, compartment-enriched RFUs, and temporal divergence.
- Technical robustness uses the frozen fractions and seeds above. Donor
  robustness uses leave-one-donor-out, 1,000 donor-level bootstrap replicates
  with seed 20260824, and exact donor-label permutation when the complete
  permutation space is feasible; otherwise 10,000 fixed-seed permutations are
  used.
- The six-volunteer cohort is a deep repeated-measures methods demonstration,
  not population-level aging evidence.

## Comparator freeze

- Every method receives identical samples, candidate sets, donor/time
  restrictions, and train/test partitions.
- Exact CDR3 and clonotype representations are count/proportion matrices;
  clonotypes are sample-namespaced if identifiers are not globally defined.
- V and J alleles are stripped before usage aggregation. CDR3 length is amino-acid
  length. Diversity summaries use the same receptor rows and declared abundance
  weights as the RFU representation.
- Edit-distance grouping uses connected components at amino-acid Levenshtein
  distance at most one and is attempted only for at most 2,000 unique sequences,
  the previously frozen quadratic-work bound. Above that bound it is reported
  unavailable, not silently replaced.
- Results are interpreted as complementary performance; RFU is not required or
  expected to win every statistic.

## VDJdb freeze

These choices are prespecified but unexecuted while `VDJDB_PATH` and
`VDJDB_RELEASE` are unset.

- The reference must be the explicit 2026-06-03 content release unless this
  document is amended before acquisition for a documented reason. The file
  SHA256, row counts, chain counts, fields, score distribution, and duplicate
  sequence-antigen records are recorded before matching.
- Analyses retain human TRB records with usable CDR3 and antigen labels. Matching
  is performed both by exact CDR3 and by exact CDR3 plus normalized V gene, with
  V alleles stripped consistently on both sides.
- Both nearest and threshold-qualified RFU policies are run on Wells 25k,
  GSE190905, and GSE157007.
- The primary evidence unit and weighting are `unique_sequence_antigen`: a
  repeated sequence-antigen record contributes once. Fractional allocation of
  ambiguous sequence labels is primary; exclusion of ambiguous sequences is a
  required sensitivity analysis.
- All validated evidence-score records are primary. Numeric score cutoffs of at
  least 1, 2, and 3 are sensitivity analyses only if the acquired release
  exposes and documents that scale; unavailable cutoffs are recorded rather
  than inferred.
- RFU summaries require at least two independently matched distinct sequences
  per RFU. Report match fraction, matched RFUs, antigen richness, purity,
  normalized entropy, same-antigen pair fraction, ambiguity, and cross-sample
  prevalence. No RFU is described as antigen-specific.
- The primary null endpoint is same-antigen pair fraction. Secondary endpoints
  include purity and normalized entropy. Size-preserving label permutations use
  1,000 permutations and seed 20260824 for unrestricted, CDR3-length-stratified,
  TRBV-stratified, and feasible TRBV-plus-length-stratified nulls.
- Grouping comparators are TRBV, CDR3 length, TRBV plus length, size-matched
  random partitions, and the frozen feasible edit-distance groups. Empty or
  singleton strata are retained according to the declared statistic rather than
  pooled after viewing results.

## Held-out freeze

GSE157007 remains governed by
[`heldout_gse157007_preregistration.json`](heldout_gse157007_preregistration.json):
the RFU reference, threshold, chain, weighting, metrics, comparator definitions,
metadata harmonization, and primary endpoints are not tuned after outcome
inspection. One sample per donor precludes within-donor retrieval in this
cohort.

## Change control

Any later change must state whether it is an implementation correction or a
scientific-definition amendment, explain why, record the date, and preserve the
previous result. Outcome inconvenience is not a valid reason to alter a primary
definition.
