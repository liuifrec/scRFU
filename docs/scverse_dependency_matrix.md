# scverse dependency compatibility matrix

Validated on 2026-08-31 in the project Python 3.10 environment.

| Dependency | Declared boundary | Tested current version | Native AIRR result |
|---|---:|---:|---|
| AnnData | `>=0.10` | 0.11.4 | PASS: X=None, H5AD, slicing, concat guard |
| MuData | optional `>=0.2.3` | 0.3.10 | PASS: modality routing, H5MU, slicing |
| Scirpy | optional `>=0.13` | 0.22.4 | PASS: current AIRR records and `index_chains` |
| Awkward | supplied by Scirpy | 2.13.0 | PASS: aligned jagged records and serialization |
| Scanpy | optional `>=1.9.3` | 1.11.5 | PASS: summary grouping/UMAP tutorial smoke |
| pandas | `>=2.0` | 2.3.3 | PASS |
| NumPy | `>=1.23` | 2.2.6 | PASS |

The coherent minimum native stack was installed in isolation and passed the
synthetic H5MU tutorial: AnnData 0.10.0, MuData 0.2.3, Scirpy 0.13.0, Awkward
2.1.0, Scanpy 1.9.3, pandas 2.0.0, and NumPy 1.23.5. The initial literal extras
bounds (`mudata>=0.2`, `scanpy>=1.9`) were tightened because Scirpy 0.13 itself
requires MuData 0.2.3 and Scanpy 1.9.3. Older DataFrame-like AnnData inputs
remain available through the 0.4 compatibility API.

No upper pins are imposed. The optional CI job resolves current compatible
versions and exercises native round trips. Future upstream AnnData/Awkward
storage changes are contained by `schema_version` and
`validate_scrfu_schema`.
