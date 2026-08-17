# Longitudinal RFU methods

The repeated-measures contract requires explicit sample, donor, and time keys.
Each sample maps to one donor and one timepoint; optional compartment,
phenotype, condition, tissue, and batch labels must also be internally
consistent. Times are numeric or explicitly ordered categoricals. Missing visits
are represented in a mask and never imputed.

`rfu_longitudinal_matrix()` reuses RFU pseudobulk semantics: one biological
sample is one row; assignment is nearest or threshold-qualified; weighting is
cell or unique sequence; normalization is count, proportion, counts per 1,000,
or CLR. Its donor/time index is a representation, not aggregation across
different samples.

Pairwise methods support cosine, Jaccard, weighted Jaccard, Bray–Curtis
dissimilarity, and Jensen–Shannon distance. Undefined zero-vector results remain
undefined. The tidy table records donor/time/compartment relations and interval;
sample pairs are not independent replicates.

Donor retrieval excludes the query and reports correct-donor rank, top-k match,
reciprocal rank, and candidate count. Parameters must be frozen outside the
evaluation cohort. Dynamics labels use declared abundance, coverage, appearance,
disappearance, fold-change, and pseudocount thresholds. They are descriptive
trajectory classes, preserve original trajectories/missingness, and do not imply
population change.

Bootstrap and permutation utilities resample or relabel entire donor blocks,
retaining their repeated observations. Percentile intervals and empirical nulls
do not turn a small cohort into population-level evidence.
