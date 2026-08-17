# Receptor Cache

Cache schema version `2` is expression-free and dataset-independent:

```text
cache/
  receptors.tsv.gz
  obs_metadata.tsv.gz
  preparation_manifest.json
```

The manifest records schema versions, adapter identity/configuration, source
format and sampled fingerprint, source dimensions when known, metadata fields,
row/cell/chain/productive counts, extraction time, scRFU/runtime versions,
adapter QC, and file SHA-256 checksums. Loading does not require the source.
Optional validation reports source state as unchanged, changed, or unavailable;
the original absolute path is not required to match.

Legacy Wells caches (`tcr_ir.tsv.gz`, schema 1) remain readable through
`scrfu.wells.load_wells_receptor_cache()`. Use
`scrfu.io.migrate_wells_receptor_cache(old, new)` or
`scrfu migrate-receptor-cache old --outdir new` for explicit migration. Generic
cache APIs never silently reinterpret legacy schemas.
