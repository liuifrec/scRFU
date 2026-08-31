# Configuration

scRFU resolves configuration in this order: an explicit Python function
argument or CLI option, then its documented environment variable when one
exists, then the documented default. There are no developer-path fallbacks.

| Setting | Explicit interface | Environment | Default |
|---|---|---|---|
| Official RFU checkout | `rfu_dir` / `--rfu-dir` | `RFU_DIR` | unconfigured; assignment raises an actionable error |
| RFU work/cache directory | `workdir` / `--workdir` | none | `.scrfu` under the current working directory |
| VDJdb reference | `path` argument in VDJdb loaders/workflows | `VDJDB_PATH` in examples that opt in | unconfigured; never downloaded silently |
| VDJdb release | explicit release argument/workflow option | `VDJDB_RELEASE` in examples that opt in | required for real-reference provenance |
| Workers | `max_workers` / `--max-workers` | none | `1` |
| Chunk size | `chunk_size` / `--chunk-size` | none | unchunked single backend call |

Native AIRR routing is always explicit in Python: `airr_mod="airr"`,
`airr_key="airr"`, `chain_idx_key="chain_indices"`, and
`key_added="scrfu"` are the current Scirpy-compatible defaults. `airr_mod` is
used only for MuData. Chain indices are optional and affect only an explicitly
requested `primary_vdj` cell summary; they never affect RFU assignments.

`RFU_DIR` must point to a user-supplied upstream RFU checkout containing the
required public RFU files. scRFU does not bundle those assets. Runtime manifests
may contain full local paths for troubleshooting; shareable evidence indexes
record filenames and hashes without machine-specific paths.

The scientific defaults—standard official RFU mode, exact-CDR3 identity, and
threshold 0.6—are frozen. Operational settings such as worker count and chunk
size must not alter assignments.
