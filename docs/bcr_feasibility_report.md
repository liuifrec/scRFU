# Experimental BCR feasibility report

Decision date: 2026-08-25. **Decision: NO-GO for a frozen BCR receptor-state
reference.** BCR canonicalization, QC, heavy/light pairing, state features, and
the structured feature matrix remain experimental utilities. No BCR
functional-unit API or reference object was constructed.

## Evidence and bounded experiment

GSE219098 was the development dataset (177,663 cells; three donors with donor
labels available for 70.0% of cells). GSE266519 was the frozen validation
dataset (16,143 cells). Construction used receptor fields only and never used
donor, phenotype, vaccine response, disease, age, timepoint, or cell-state
labels.

The acquired tables lack mutation frequency, germline identity, and inferred
clonal-family identifiers. Therefore maturation-aware and clonal-family-aware
candidates could not be specified without fabricating data. Three exact,
weight-free baselines were evaluated instead:

| Outcome-independent candidate | Development groups | Singleton groups | 50% holdout coverage | Frozen GSE266519 coverage | Development donor purity |
|---|---:|---:|---:|---:|---:|
| Heavy CDR3 exact | 160,263 | 95.62% | 8.84% | 0.236% | 99.83% |
| Heavy V/J + CDR3 exact | 160,654 | 95.75% | 8.60% | 0.0128% | 99.92% |
| Paired heavy/light exact | 160,845 | 96.47% | 6.69% | 0% | 99.99% |

The exact representations are legitimate comparators, but they are not novel
receptor-state groups. Their near-total singleton structure, poor half-sample
coverage, negligible frozen cross-dataset coverage, and extreme donor
concentration do not support transferability. Introducing a fuzzy distance
would require unvalidated edit costs or feature weights, which this sprint did
not tune against validation labels.

## Prespecified gate

| Criterion | Result | Evidence |
|---|---|---|
| Deterministic construction | PASS for exact baselines | Stable exact tuple keys and five fixed split seeds |
| Clear receptor-specific feature definition | PARTIAL | Exact heavy/paired definitions are clear; no justified non-exact distance |
| Reasonable held-out reference coverage | FAIL | 0–0.236% across candidates |
| Not overwhelmingly donor-specific | FAIL | Weighted donor purity 99.83–99.99% among donor-labelled development cells |
| Not merely exact clonal families | FAIL / untestable | Exact keys are clonotype-like; inferred family identifiers are absent |
| Adequate downsampling stability | FAIL | 50% holdout coverage 6.69–8.84% |
| Cross-dataset mapping without refitting | FAIL | At most 37 validation cells covered by the heavy-CDR3 baseline |
| Interpretable relationship without fitting labels | PARTIAL | Exact groups are highly isotype/clonotype concentrated, but that is not functional convergence |
| Documented missing-data behavior | PASS | Required-field missingness gives an explicit unassigned state |
| Non-trivial relative to one baseline variable | FAIL | Only exact baseline tuples were defensible with available fields |

## Recommendation

Retain the experimental preprocessing and feature-extraction surface. Do not
expose BCR functional units, do not reuse TCR centroids, and do not claim
functional convergence. A future gate would require a public development and
independent cohort with compatible germline/SHM and family annotations, a
prespecified biologically justified distance, donor-aware construction, and
frozen cross-cohort validation. This future work is optional and is not a TCR
scRFU release blocker.

The authoritative external report is
`bcr_representation_feasibility.json` (SHA256
`fa692f99748df4bf3c6949bcc449342bce1bbbef5a9ef5af87d926574389407d`).
The path is intentionally relative; public data and derived tables remain
outside the repository.

