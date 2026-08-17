# Memory-Efficient Wells Receptor Preparation

> Compatibility note: this page describes legacy Wells cache schema 1.
> Dataset-independent cache schema 2 uses `receptors.tsv.gz`; see
> [receptor_cache.md](receptor_cache.md). Legacy reading remains supported and
> migration is explicit and non-overwriting.

## Why ordinary H5AD loading is expensive

The Wells atlas H5AD is compressed on disk, but ordinary
`anndata.read_h5ad()` materializes its arrays. In the inspected 610,429-cell,
35,475-feature atlas, the CSR arrays for `X` occupy about 11.19 GiB after
decompression and `raw/X` occupies about 11.28 GiB. All `uns` entries add about
1.28 GiB of logical array data. The `uns["TCR_IR"]` dataframe itself accounts
for only about 0.064 GiB.

Slicing after ordinary loading cannot avoid that cost: both expression matrices
and all `uns` values already exist in memory before the first 1,000 observations
are selected. Backed AnnData avoids materializing `X` and `raw/X`, but AnnData
still reads `uns` eagerly. Backed mode therefore improves memory substantially
but still reads unrelated BCR and CITE-seq data that receptor preparation does
not need.

## Targeted reader

`scrfu.wells.read_wells_receptors_h5ad()` opens the HDF5 container directly and
reads only:

- the observation index;
- explicitly selected `obs` columns; and
- the row-aligned `TCR_IR` dataframe from `uns` or, as a compatibility fallback,
  `obsm`.

It does not read `X`, `raw`, `layers`, `obsp`, or unrelated `uns` values. It does
not assume that an arbitrary group is a dataframe. At runtime it validates the
AnnData root encoding, encoded-dataframe type, declared index, column order,
column presence, categorical representation, row counts, and exact `obs` to
`TCR_IR` index alignment. Unsupported encodings raise
`UnsupportedWellsH5ADLayout` with the failing logical element.

## Compact cache

Prepare an expression-free cache with:

```bash
scrfu prepare-wells /path/to/wells_atlas.h5ad \
  --output-dir /path/to/wells_receptor_cache \
  --obs-column donor_id \
  --obs-column library_id \
  --obs-column cell_type
```

For a restricted smoke cache, add `--max-cells 1000`. The cache contains:

```text
wells_receptor_cache/
  tcr_ir.tsv.gz
  obs_metadata.tsv.gz
  preparation_manifest.json
```

The manifest records cache schema version, scRFU version, UTC extraction time,
source atlas dimensions, selected cell and TCR row counts, selected metadata
columns, TCR container, and source identity. Source identity includes the
resolved source path, byte size, nanosecond modification time, and a SHA-256
digest over the file size plus fixed first, middle, and last 1 MiB samples. This
lightweight fingerprint avoids hashing an entire multi-gigabyte atlas while
detecting metadata changes and sampled content changes. Supplying the source
again when loading a cache enforces fingerprint equality:

```bash
python examples/wells_atlas_workflow.py \
  --input /path/to/wells_receptor_cache \
  --source-h5ad /path/to/wells_atlas.h5ad \
  --outdir /path/to/results \
  --rfu-dir /path/to/RFU
```

Each cache table also has a SHA-256 checksum in the manifest. Cache loading
requires exact cell-index alignment and matching recorded row counts. The cache
never contains an expression matrix.

## Limitations

- The targeted reader supports AnnData encoded dataframes with one-dimensional
  datasets, categoricals, and standard nullable scalar encodings. It fails
  explicitly for a new or unsupported encoding instead of guessing.
- A lightweight source fingerprint is not a full-file cryptographic hash. Use a
  separately managed full SHA-256 when byte-for-byte archival identity is
  required.
- `--write-annotated` requires original H5AD input and is an explicit opt-in to
  expression loading. It is unavailable for receptor-cache input because the
  cache deliberately contains no expression data.
- A prefix selected with `--max-cells` is a smoke-test cohort, not a statistically
  representative atlas sample unless the source ordering makes it so.
