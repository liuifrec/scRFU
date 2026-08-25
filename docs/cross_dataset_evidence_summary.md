# Cross-dataset evidence summary

Updated on 2026-08-25. Numbers below are audited in
[`manuscript_claim_audit.md`](manuscript_claim_audit.md). Dataset roles are not
interchangeable: Wells is a public single-cell development/stress demonstration,
GSE190905 is an independent public paired technical validation, and GSE157007
is preregistered held-out validation. The independent software paper does not
depend on a private cohort.

| Dataset | Fixed role | Sample/donor design | Receptor rows | Unique CDR3s | Threshold coverage | RFU richness | Threshold-qualified RFUs | Robustness evidence | Comparator evidence | Antigen evidence | Longitudinal evidence |
|---|---|---|---:|---:|---:|---:|---:|---|---|---|---|
| Wells full receptor-only atlas | Development demonstration and public paired single-cell stress test | 610,429 source cells; 303,088 productive primary-TRB receptor rows; 209 libraries; 21 receptor-bearing donors | 303,088 | 192,675 | 77.18% | 4,996 nearest RFUs | 4,928 | Exact fresh/resume reconstruction; full downstream library/donor passes; prior bounded cell/sequence/abundance robustness and chunk/order sensitivity | Exact CDR3, V, CDR3 length, Shannon and Simpson; full library pseudobulk/coupling | CDR3: 4,701 matched sequences and 2,084 matched nearest RFUs; CDR3+V: 934 RFU sequences represented by 985 strict query variants and 716 matched nearest RFUs | None; cross-sectional atlas |
| GSE190905 | Independent frozen-reference and public paired technical validation | 12 paired pre/post samples from 6 treated patients | 27,655 | 16,391 | 78.31% | 4,465 nearest RFUs | 3,940 | Within-person versus between-person descriptive similarities and leave-one-timepoint-out retrieval; no outcome testing | RFU, exact CDR3, clonotype, V, J, length, Shannon and Simpson with identical candidates; nearest-RFU cosine top-1/top-3/MRR 0.667/0.833/0.764 | CDR3: 616 matched sequences and 506 nearest RFUs; strict CDR3+V: 93 sequences and 91 nearest RFUs; strict coherence is sparse | Public paired methods demonstration; nearest RFU mean cosine 0.491 within versus 0.063 between donors; no population inference |
| GSE157007 | Preregistered held-out aging/frailty transfer validation | 17 cross-sectional samples from 17 donors in 4 age groups | 60,125 | 43,082 | 76.88% | 4,898 nearest RFUs | 4,618 | Three-seed deterministic sample-vector subsampling; nearest-RFU mean cosine 0.936 at 50% and 0.976 at 75% | RFU, exact CDR3, namespaced sample-local clonotype, V, J, length, Shannon and Simpson; edit-distance skipped by preregistered scale limit | CDR3: 1,635 matched sequences and 1,074 nearest RFUs; strict CDR3+V: 294 sequences represented by 302 query variants and 264 nearest RFUs | Undefined by design: one sample per donor |

## Interpretation boundary

The three public datasets demonstrate deterministic execution, full-atlas
receptor/downstream scale, frozen-reference coverage, transfer across distinct
cohort designs, public paired donor structure and linkage to single-cell
phenotype annotations. The VDJdb analysis adds external antigen-label coherence
evidence but does not establish antigen specificity. No result establishes a
population-level aging effect or universal comparator superiority.
