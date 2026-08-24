# Held-out transfer validation plan

A held-out cohort is declared before its RFU outcomes are inspected. The
machine-readable manifest records development, validation, and held-out cohort
roles; frozen reference identifier; thresholds and eligibility; input hashes;
evaluation metrics; and code version.

The held-out cohort cannot be used for reference fitting, threshold selection,
metric definition, robustness choices, comparator tuning, phenotype category
selection, or figure cutoffs. Metadata mappings are explicit and archived.

Primary candidate metrics are reference coverage, RFU richness/prevalence,
sample-vector stability, donor retrieval where repeated measures exist, and
predeclared comparator performance. Secondary exploratory outputs must be
labelled as such. Failure due to missing fields or low coverage is itself
reported, not repaired by target-cohort refitting.

Acceptance requires a manifest created before evaluation, verified hashes, no
role overlap, a frozen reference match, source tables, all undefined outcomes,
and a report of cohort shift and missingness.

GSE157007 was registered for the held-out role in
[`heldout_gse157007_preregistration.json`](heldout_gse157007_preregistration.json).
Its 17 receptor-file hashes were frozen before RFU assignment. The cohort is
cross-sectional with one sample per donor, so donor retrieval is undefined and
is not replaced by an age-group classifier. Sample-level subsampling stability,
reference coverage, RFU richness, and the prespecified comparators are the
applicable endpoints.
