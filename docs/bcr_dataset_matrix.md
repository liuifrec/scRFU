# Public BCR dataset matrix

Verified 2026-08-25 from the linked GEO, 10x Genomics, and OAS records. Blank or
“not verified” entries are deliberately not inferred. Public labels are
evaluation metadata only; they must not be used to construct a receptor-state
reference.

| Dataset | Design | Receptors and fields | Biological metadata | Access / sprint use |
|---|---|---|---|---|
| [GSE219098](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE219098) | Three vaccinees, repeated sampling through six months; 81 GEX/ADT/VDJ records | Single-cell paired IGH/IGK/IGL; V/D/J/C, CDR3aa/CDR3nt and 10x clonotype available. In the acquired tables, isotype is 49.16% complete; SHM, germline identity, and inferred clonal family are absent. | Donor demultiplexing, time, sorted spike-binding/non-binding and plasmablast labels | Public GEO/SRA. Acquired development dataset; processed RAW archive 165,007,360 bytes, SHA256 `32d1ffa1e6a3d73f6082aed6a09ed78eca28b8f2b51d7daa2d538709710b4eb7`. |
| [GSE244297](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE244297) | Five volunteers before and weeks 1–3 after vaccination; one excluded by the source study; 24 records | Eight single-cell 10x VDJ libraries; exact local field completeness not measured | Spike-reactive activated B-cell/plasmablast and vaccination-time labels | Public GEO/SRA; not acquired in this bounded sprint. Processed archive advertised as about 117.7 MB. |
| [GSE266519](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE266519) | Four transplant recipients and three healthy controls, with pre/post vaccination and some post-infection samples | Single-cell paired IGH/IGK/IGL; V/D/J/C, CDR3aa/CDR3nt and 10x clonotype. Isotype 48.84% complete; SHM, germline identity, and inferred family absent. | Condition and sampling labels exist in the study, but compatible per-cell donor metadata were not exported with the aggregate table | Public GEO/SRA. Acquired validation dataset; contig gzip SHA256 `ce615c2ca41ed460617328f1b0280d8ad348974e6636eb156e6174cd7a562f3c`. |
| [10x healthy donor pre/post influenza vaccination](https://www.10xgenomics.com/cn/datasets/human-b-cells-from-a-healthy-donor-pre-and-post-flu-vaccination-v-2-2-standard-5-0-0) | One healthy donor, pre-vaccination and day 14; 18,115 cells | Cell Ranger paired BCR with V/D/J/C, CDR3 and clonotype; SHM/family not advertised | Timepoint | CC BY 4.0; not acquired. Useful for adapter examples, not donor transfer. |
| [OAS paired collection](https://opig.stats.ox.ac.uk/webapps/oas/documentation_paired/) | Heterogeneous contributed studies | Paired annotated heavy/light CSVs; fields vary by study and can include germline alignments | Subject, disease, vaccination, tissue, cell type, and longitudinal metadata where supplied | OAS states CC BY 4.0. No collection was downloaded; a future use must freeze a bounded source study and hashes first. |

## Acquired-data QC

| Metric | GSE219098 | GSE266519 |
|---|---:|---:|
| Input contigs | 800,220 | 34,477 |
| Productive contigs | 379,741 | 34,477 |
| Cells | 177,663 | 16,143 |
| IGH / IGK / IGL | 187,133 / 117,043 / 75,565 | 16,852 / 10,636 / 6,989 |
| Paired / heavy-only / light-only | 171,642 / 3,087 / 2,934 | 15,494 / 177 / 472 |
| V / J completeness | 100% / 100% | 100% / 100% |
| Constant-region completeness | 99.49% | 99.94% |
| 10x clonotype completeness | 100.00% | 100.00% |
| Mutation / germline identity / inferred-family completeness | 0% / 0% / 0% | 0% / 0% / 0% |
| Exact-CDR3 duplicate row fraction | 43.20% | 35.58% |
| Malformed or missing CDR3aa | 0 | 0 |

Both runs retained all productive contigs, selected primary heavy and light
chains deterministically using productivity, UMI/read counts, and stable source
order, and preserved source IDs. The feature output has exactly one row per
cell. Multiple valid candidates remain in the canonical contig table and are
reported by QC rather than silently discarded.

