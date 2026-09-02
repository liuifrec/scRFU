# scverse representation benchmark

This benchmark compares representations on the same 12 GSE190905
patient-state samples (six repeated donors) without fitting or tuning on donor
identity. Genuine Scirpy 0.22.4 clonotypes were generated with nucleotide
identity distance, all receptor arms, primary-only dual-IR handling, and the
documented default within-receptor-type definition. All 27,655 receptor-bearing
cells entered every count representation.

## Representation properties

| Representation | Features | Matrix sparsity | Effective dimension |
|---|---:|---:|---:|
| RFU | 4,465 | 0.749 | 1,582.2 |
| Exact CDR3 | 16,391 | 0.913 | 4,619.4 |
| Scirpy clonotype | 16,937 | 0.913 | 4,936.9 |
| TRBV + TRBJ | 561 | 0.391 | 232.0 |
| CDR3 length | 15 | 0.178 | 6.6 |

## Repeated-donor results

| Representation | Top-1 | Top-3 | MRR | Within-donor cosine | Between-donor cosine | 50% subsample cosine | 75% subsample cosine |
|---|---:|---:|---:|---:|---:|---:|---:|
| RFU | 0.667 | 0.833 | 0.764 | 0.491 | 0.063 | 0.946 | 0.979 |
| Exact CDR3 | 0.667 | 0.750 | 0.771 | 0.477 | 0.019 | 0.939 | 0.977 |
| Scirpy clonotype | 0.667 | 0.667 | 0.733 | 0.476 | 0.019 | 0.934 | 0.976 |
| TRBV + TRBJ | 0.667 | 0.833 | 0.778 | 0.591 | 0.334 | 0.975 | 0.991 |
| CDR3 length | 0.417 | 0.833 | 0.653 | 0.941 | 0.918 | 0.999 | 1.000 |
| Diversity summaries | 0.000 | 0.417 | 0.297 | not a count-space comparison | not a count-space comparison | not evaluated as a count vector | not evaluated |

The result is complementary, not a universal RFU win. V/J was the best of
these representations by MRR and downsampling cosine; exact CDR3 slightly
exceeded RFU MRR. RFU used substantially fewer features than exact CDR3 or
Scirpy clonotypes, had slightly higher subsampling stability than those two,
and preserved the same top-1 retrieval. Very low-dimensional summaries are
stable partly because they discard receptor detail.

Audited external sources are `representation_summary.tsv`,
`retrieval_summary.tsv`, `pairwise_summary.tsv`, and
`downsampling_summary.tsv` under the GSE190905 Scirpy-comparator evidence root;
their hashes are recorded in `docs/source_data_inventory.md`.
