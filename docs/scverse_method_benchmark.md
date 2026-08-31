# scverse representation benchmark organization

This document reorganizes existing validated outputs; it does not rerun or
retune scientific analyses. Real tables remain external and are indexed by
`docs/evidence_index.md`.

| Representation | Dimensionality/sparsity | Downsampling/stability | Donor/longitudinal | Frozen transfer | Runtime/memory |
|---|---|---|---|---|---|
| RFU | Sample × frozen RFU, lower dimensional than exact CDR3 | Wells robustness and GSE157007 50%/75% cosine evidence | GSE190905 donor retrieval and within/between similarity | GSE190905 and preregistered GSE157007 coverage | Full Wells and bounded benchmarks |
| Exact CDR3 | Very sparse cohort-specific sequence space | Existing comparator/downsampling tables | Existing GSE190905 comparator | Exact transfer limited by sequence recurrence | Comparator source tables |
| Scirpy clonotype | Dataset-defined clonotype identity where available | Same sample restrictions as comparator runs | Public repeated-donor comparator | Not a frozen external reference | Existing comparator source tables |
| V/J usage | Low-dimensional conventional baseline | Existing robustness tables | Existing donor comparator | Directly harmonizable after gene normalization | Existing comparator source tables |
| CDR3 length/diversity | Very low-dimensional summaries | Existing robustness tables | Existing comparator | Directly harmonizable | Existing comparator source tables |
| Sequence-distance groups | Data-dependent clusters | Existing deterministic edit-distance comparator | Existing comparator | Requires fixed grouping definition for transfer | Existing comparator source tables |

The benchmark does not assert universal RFU superiority. The supported method
claim is complementary: RFUs provide an interpretable frozen-reference state
space with measurable coverage, stable bounded subsampling behavior, and
cross-cohort mapping, while conventional representations retain advantages for
other tasks.

Exact values and hashes must be taken from the external manifests, not this
organizational summary.

