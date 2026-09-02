# Public longitudinal TCR candidates for scverse validation

Audit date: 2026-09-01. This is a public-data methods-validation inventory, not
a biological outcome analysis. No private cohort was searched or used.

## Ranked candidates

| Rank | Accession | Repeated design | Processed receptor availability | Decision |
|---|---|---|---|---|
| 1 | [GSE345124](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE345124) / PRJNA1518857 | 18 donors; days 0, 42 and 63; 13 donors have all three visits and five have two; paired CD4/CD8 tables at each observed visit | 98 public processed Adaptive TCR TSV.GZ files; 164,392,960-byte archive; no raw sequence release | Strongest verified deeper repeated-measures candidate. Acquisition/QC complete; full RFU/Scirpy execution deferred because it requires 2,258,994 distinct eligible CDR3 queries and 3,038,924 in-frame rows. No convenience subsampling was substituted. |
| 2 | [GSE190905](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE190905) / PRJNA788905 | Six donors with paired pre/post samples | Compact processed single-cell TCR and metadata tables | Completed bounded longitudinal demonstration, including genuine Scirpy clonotypes, within/between similarity, retrieval, and downsampling. Two visits limit dynamics claims. |
| 3 | PRJNA602091 | Repeated donors with two visits separated by about nine years | Raw sequencing rather than a compact processed receptor table | Scientifically relevant but does not meet the preferred three-timepoint/processed-table criteria for this sprint. |

## GSE345124 verified QC

The downloaded archive SHA256 is
`5e7df6742d41da5a7cd71c1b1ffbc60f3a368bc64ef1b783c77c56be049bf1a9`.
Its 98 tables use two documented Adaptive export schemas (80 legacy and 18
current). They contain 3,740,032 rows, of which 3,038,924 are in frame and
3,032,688 meet the canonical CDR3 shape gate. All 49 observed donor-timepoint
combinations have both CD4 and CD8 files. The external QC manifest is indexed
in `docs/evidence_index.json`; no receptor rows are stored in Git.

The dataset title and sample descriptions contain treatment/disease context.
Any later execution is a representation-method validation and must not turn
those labels into fitted features or make unprespecified clinical claims. A
publication identifier was not present in the retrieved GEO family metadata at
the audit date, so it remains an explicit citation caveat rather than an
invented reference.
