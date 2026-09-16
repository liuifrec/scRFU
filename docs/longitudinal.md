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

## Change hidden by fixed grouping

`scrfu.tl.multiscale_repertoire_change(before, after, grouping)` accepts two
nonnegative receptor-mass Series and one fixed receptor-to-group Series. It
returns normalized receptor and group changes, total-variation distances and
`aggregation_cancellation = d_clone - d_group`. The group distance cannot
exceed the receptor distance. This standard contraction identity is not a
measure of preserved antigen recognition or functional recovery.

Unmapped receptors raise by default. With `unmapped="condition"`, both distances
are conditional on mapped receptors; original and retained masses/coverage
remain in the summary. Empty mapped samples yield unknown distances, including
when both are empty. Missing visits must be handled in the design table, not
passed as observed empty samples. Counts are never reconstructed from normalized
values. Callers retain responsibility for clone identity and sampling units.

`scrfu.tl.permute_fixed_groups(grouping, strata=features, random_state=...)`
returns a fixed shuffled map preserving group-by-stratum counts over the
supplied receptor universe. Reuse the map at every visit. Report strata that
cannot exchange labels and do not call this a calibrated biological null.
The [radiotherapy application specification](../manuscript/radiation_methods_specification.md)
defines one use with explicit donor-level reporting and empirical cell
subsampling, without changing assignment.
