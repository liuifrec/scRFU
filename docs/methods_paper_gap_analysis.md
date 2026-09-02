# Methods-paper gap analysis

Audit date: 2026-09-01. This document evaluates evidence; it is not manuscript
prose. “Demonstrated” means a completed analysis with an externally hashed
source table, not merely an implemented function.

| Potential method claim | Status | Exact evidence | Analysis unit and caveat |
|---|---|---|---|
| 1. Exact fidelity to canonical RFU assignment | DEMONSTRATED | Official-RFU parity manifest indexed in `evidence_index.json`; 10/10 eligible rows, zero ID/label/threshold/order mismatches and maximum score difference 0 | Adversarial receptor row; tested official reference/artifacts only |
| 2. Scalable execution on large single-cell repertoires | DEMONSTRATED | Full Wells: 610,429 source cells, 192,675 unique queries, 418 s backend work, 3.42 GB peak RSS, 1.33 s backend resume; native 25k/100k/250k source tables | Technical runtime on the Linux development host; not a hardware-independent speed claim |
| 3. Native AnnData/MuData/AIRR interoperability | DEMONSTRATED | Native schema tests plus public Scirpy `wu2020_3k` manifest: 7,544 chains, 2,931 eligible TRBs, zero table/native mismatches | AIRR chain / observation; official RFU execution validated on Linux |
| 4. Stable serialization, subsetting and guarded concatenation | DEMONSTRATED | H5AD/H5MU round-trip, slicing, compatible concat, incompatible reference/schema rejection tests; bounded native reload mismatch count 0 | Technical object behavior; Awkward-in-AnnData support remains upstream experimental |
| 5. Compression relative to exact receptor identity | DEMONSTRATED | `representation_compression.tsv`: Wells 192,675 CDR3s to 4,996 RFUs; GSE190905 16,391 to 4,465; GSE157007 43,082 to 4,898; wu2020_3k 2,295 to 1,759 | Dataset-level representation property; compression is not biological superiority |
| 6. Frozen-reference transfer across independent datasets | DEMONSTRATED | GSE190905 independent transfer, preregistered GSE157007 held-out transfer, and `cross_dataset_feature_sharing.tsv` | Cohort/sample; no RFU refit or held-out tuning; differing cohort purposes limit biological pooling |
| 7. Improved or complementary longitudinal donor representation | PARTIAL | GSE190905 genuine comparator: RFU top-1/top-3/MRR 0.667/0.833/0.764; within/between cosine 0.491/0.063 | Twelve samples from six donors and two visits; supports complementarity but not deep trajectory dynamics or universal improvement |
| 8. Single-cell phenotype integration | DEMONSTRATED descriptively | Full/bounded Wells phenotype-coupling and pseudobulk source tables; native synthetic MuData tutorial | Cell phenotype and library/donor summaries; no cell-level inferential p-values or causal claim |
| 9. External antigen-annotation coherence | DEMONSTRATED narrowly | Pinned VDJdb 2026-06-03 24-combination analysis and nulls; Wells nearest/fractional CDR3 observed 0.1128 versus null 0.0646, empirical p 0.000999 | Distinct matched sequence annotation; external coherence only, never antigen specificity |
| 10. Complementarity to Scirpy clonotype analysis | DEMONSTRATED | Genuine GSE190905 Scirpy clonotypes and identical candidate/downsampling sets in `scverse_method_benchmark.md` | Sample representation; RFU does not win every endpoint |

## Unsupported claims excluded

- RFU universally outperforms exact CDR3, Scirpy clonotypes, or V/J summaries.
- Individual RFUs are antigen-specific.
- Two-timepoint GSE190905 establishes deep temporal RFU dynamics.
- Cross-cohort feature reuse proves biological equivalence between cohorts.
- Official RFU/R execution is supported on macOS or Windows; only Python-only
  functionality is CI-tested there.

## Smallest remaining scientific experiment

The only material evidence gap for a stronger longitudinal claim is a complete,
prespecified run on a public cohort with at least three visits. GSE345124 is now
verified as the strongest candidate (18 donors, 13 with three visits), but its
2,258,994 unique eligible CDR3 queries place it outside this bounded sprint.
This does not block the narrower transferability and two-visit repeated-donor
claims. It does block any claim about persistent/expanding/contracting RFU
trajectories across three or more visits.
