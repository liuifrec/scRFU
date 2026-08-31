# scverse readiness scorecard

Assessment date: 2026-08-31. Detailed machine-readable evidence is in
`docs/scverse_readiness_scorecard.json`.

| Domain | Status | Basis |
|---|---|---|
| Native AnnData / MuData | PASS | Current AIRR Awkward and modality tests |
| AIRR semantics / chain alignment | PASS | Every source chain retained and aligned |
| Scirpy interoperability | PASS | Current `index_chains` and end-to-end flow |
| Serialization / subsetting | PASS | H5AD and H5MU round trips and slicing |
| Concatenation | PASS | Compatible references preserved; incompatible references rejected |
| X independence / memory safety | PASS | Receptor-only code path; targeted large-H5AD path retained |
| API / schema / provenance | PASS | `assign_rfu`, schema 1.0, portable hashes, validator |
| Documentation / tutorial | PASS | Sphinx sources and synthetic MuData tutorial |
| Coverage | PARTIAL | Measured; user-visible native path is covered, legacy gaps remain |
| CI | PARTIAL | Optional scverse job configured; remote run requires user push |
| Dependency compatibility | PASS | Current and coherent minimum stacks pass the native tutorial |
| Standard registry | BLOCKED | No PyPI/conda release yet |
| Hosted documentation | PARTIAL | ReadTheDocs config ready; project not provisioned |
| Ecosystem registry | BLOCKED | Registry/docs publication and maintainer consent pending |
| Method validation | PASS | Frozen 0.4 public evidence remains applicable |

These classifications deliberately separate implementation readiness from
external service actions. The lack of a standard-registry release is a real
mandatory ecosystem blocker; it is not a scientific-method defect.
