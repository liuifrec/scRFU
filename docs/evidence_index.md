# External evidence index

This index records hashes of compact evidence manifests stored beside external
runtime outputs. It does not redistribute RFU assets, VDJdb, validation data,
or derived result tables. Paths in the JSON index are portable filenames or
relative labels; machine-specific runtime paths remain only in external
manifests.

The index is generated after external manifests are sealed. The historical
official-RFU parity result is retained with its previously audited hash if its
temporary runtime directory is no longer available; that status is stated
explicitly rather than reconstructing evidence.

| Dataset | Analysis | Manifest filename | SHA256 | Status |
|---|---|---|---|---|
| Official RFU parity fixture | canonical parity | `original_rfu_parity/summary.json` | `4adc57d3…43e5` | previously verified; temporary source unavailable for resealing |
| Wells atlas | full RFU assignment | `evidence_manifest_full_rfu.json` | `9e93719e…3ec8` | sealed |
| Wells atlas | full receptor downstream | `evidence_manifest_downstream.json` | `5411ec56…d7be` | sealed |
| GSE190905 | frozen RFU transfer | `evidence_manifest_rfu.json` | `e7bb244d…4d13` | sealed |
| GSE157007 | preregistered held-out transfer | `evidence_manifest_rfu.json` | `260abade…6eab` | sealed |
| VDJdb 2026-06-03 | external antigen-evidence validation | `evidence_manifest.json` | `c32da315…f242` | sealed |
| Synthetic scaling | pure-Python benchmark | `evidence_manifest.json` | `0025641d…07c6` | sealed |
| Wells atlas | completed-run validation | `full_run_validation_report.json` | `e613840a…c043` | sealed |
| VDJdb representative subsets | deterministic reproducibility check | `reproducibility_check.json` | `bdcd7785…fbd4` | sealed |
| GSE219098 | experimental BCR adapter and feature QC | `GSE219098/qc/run_manifest.json` | `bf7ff9b6…653b` | sealed |
| GSE266519 | experimental BCR adapter and feature QC | `GSE266519/qc/run_manifest.json` | `0044c846…e8db` | sealed |
| GSE219098 → GSE266519 | experimental BCR feasibility | `bcr_representation_feasibility.json` | `fa692f99…407d` | sealed |
| Wells atlas | native AIRR 25k scale/parity | `native_scale/25k/run_manifest.json` | `fe65d608…64c7` | sealed |
| Wells atlas | native AIRR 100k scale/parity | `native_scale/100k/run_manifest.json` | `7ae650ae…bf50` | sealed |
| Wells atlas | native AIRR 250k scale/parity | `native_scale/250k/run_manifest.json` | `e79c555b…9395` | sealed |
| GSE190905 | genuine Scirpy clonotype comparator | `GSE190905/scirpy_comparator/run_manifest.json` | `d75ee0e8…52d1` | sealed |
| Wells/GSE190905/GSE157007/wu2020_3k | representation compression and sharing | `representation_consistency/run_manifest.json` | `942ad5b0…18c0` | sealed |
| Scirpy wu2020_3k | official-RFU native H5MU interoperability | `wu2020_3k/official_run/run_manifest.json` | `3bcab11e…1cd2` | sealed |
| Scirpy wu2020_3k + VDJdb 2026-06-03 | native chain-level VDJdb linkage | `wu2020_3k/native_vdjdb_linkage/run_manifest.json` | `3e84401e…79d3` | sealed |
| GSE345124 | public longitudinal candidate acquisition/QC | `GSE345124/qc/run_manifest.json` | `3fbd7d99…0b78` | verified; frozen-reference execution deferred by bounded gate |

Each sealed external manifest includes full SHA256 values, byte sizes,
parameters, software/environment metadata, and the paths used at runtime. The
runtime paths are intentionally not reproduced in this shareable index.
