# Frozen-reference cross-cohort transfer

`FrozenRFUReference` records SHA256 values for `RFU.R`, the embedding artifact,
and centroid atlas; threshold; eligibility rule; assignment mode; receptor
chain/model; reference label; package version; timestamp; and schema version.
The immutable identifier hashes scientific reference parameters, not runtime
paths or timestamps.

`transfer_cohort()` accepts already assigned target-cohort rows and an immutable
reference. It validates an observed reference identifier when supplied, then
returns sample/optional-group coverage, score quantiles, RFU
richness/prevalence, sample pseudobulk, source rows, and provenance. It cannot
fit a reference, change a threshold, or perform expression-level batch
integration. A target result without an observed reference ID is clearly marked
unverified.

`harmonize_cohort_metadata()` maps canonical fields to source fields and values
only through caller-provided dictionaries. It reports missing and unmapped
values plus source-only fields; strict mode rejects any unmapped value. It does
not infer an ontology.

Held-out manifests freeze cohort roles, parameters, input hashes, metrics, code
version, and reference identity before evaluation. Development, validation, and
held-out roles cannot overlap.
