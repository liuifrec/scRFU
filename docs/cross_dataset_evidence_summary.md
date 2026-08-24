# Cross-dataset evidence summary

Frozen on 2026-08-24. Numbers below are audited in
[`manuscript_claim_audit.md`](manuscript_claim_audit.md). Dataset roles are not
interchangeable: Wells is a bounded development demonstration, GSE190905 is an
independent paired technical validation, GSE157007 is preregistered held-out
validation, and the governed six-volunteer cohort is the intended deep
longitudinal demonstration.

| Dataset | Fixed role | Sample/donor design | Receptor rows | Unique CDR3s | Threshold coverage | RFU richness | Threshold-qualified RFUs | Robustness evidence | Comparator evidence | Antigen evidence | Longitudinal evidence |
|---|---|---|---:|---:|---:|---:|---:|---|---|---|---|
| Wells bounded 25k | Development demonstration and public paired single-cell atlas | 25,000 sampled source cells; 12,415 receptor-bearing cells; 209 libraries; 21 receptor-bearing donors; 16 receptor-bearing cell types | 12,415 | 10,038 | 76.94% | 3,963 nearest RFUs | 3,381 | Cell/sequence subsampling, multinomial abundance resampling, threshold/chunk/order/cache sensitivity; nearest cell-subsampling median cosine 0.526/0.730/0.885 at 25/50/75% | Exact CDR3, V, CDR3 length, Shannon, Simpson; J and clonotype unavailable; edit-distance skipped by frozen scale limit | Blocked: VDJdb unset | None; cross-sectional atlas |
| GSE190905 | Independent frozen-reference technical validation | 12 paired pre/post samples from 6 treated patients | 27,655 | 16,391 | 78.31% | 4,465 nearest RFUs | 3,940 | Within-person versus between-person descriptive similarities and leave-one-timepoint-out retrieval; no outcome testing | RFU, exact CDR3, clonotype, V, J, length, Shannon, Simpson with identical candidates; nearest-RFU cosine top-1/top-3/MRR 0.667/0.833/0.764 | Blocked: VDJdb unset | Paired technical structure only; not the governed longitudinal cohort and not population evidence |
| GSE157007 | Preregistered held-out aging/frailty transfer validation | 17 cross-sectional samples from 17 donors in 4 age groups | 60,125 | 43,082 | 76.88% | 4,898 nearest RFUs | 4,618 | Three-seed deterministic sample-vector subsampling; nearest-RFU mean cosine 0.936 at 50% and 0.976 at 75% | RFU, exact CDR3, namespaced sample-local clonotype, V, J, length, Shannon, Simpson; edit-distance skipped by preregistered scale limit | Blocked: VDJdb unset | Undefined by design: one sample per donor |
| Governed six-volunteer cohort | Deep longitudinal methodological demonstration | Required design is donor×time×optional compartment; actual counts unavailable | Blocked | Blocked | Blocked | Blocked | Blocked | Prespecified cell/sequence/abundance resampling only | Prespecified identical-candidate RFU and conventional representations | Would use the same external VDJdb analysis only if configured | Blocked: explicit receptor, metadata, and output paths unset |

## Interpretation boundary

The three completed public datasets demonstrate deterministic execution,
bounded stability, frozen-reference coverage, transfer across distinct cohort
designs, and linkage to single-cell phenotype annotations. They do not establish
population-level aging effects, universal comparator superiority, antigen
specificity, or governed longitudinal dynamics. The first two missing evidence
classes require real VDJdb and explicit longitudinal inputs respectively; their
software implementations alone are not evidence.
