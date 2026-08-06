# API Contract

This document describes the intended stable public API for manuscript reviewers
and future users. scRFU is currently pre-alpha, but these functions are the
planned compatibility surface.

## Canonical receptor contract

`scrfu.pp.validate_receptor_table()` validates schema-versioned canonical
receptor rows. `scrfu.adapters.prepare_receptors()` resolves an explicit named
adapter and returns receptor rows, separate cell metadata, QC, and provenance.
Built-ins are `wells_tcr_ir`, `scirpy_airr`, `generic_airr_dataframe`, and
`cellranger_vdj`.

`scrfu.io.read_h5ad_dataframe()` and `read_h5ad_obs()` selectively read encoded
dataframes without materializing expression payloads. Generic cache schema 2 is
implemented by `write_receptor_cache()`, `read_receptor_cache()`, and
`validate_receptor_cache()`.

### `scrfu.tl.call_rfu_table`

This is the dataset-independent RFU center. It accepts canonical rows, selects
the requested chain, applies the validated eligibility/deduplication/chunking
backend, and returns `RFUTableResult(per_sequence, per_row, mapping,
provenance)`. Canonical row IDs and source provenance are preserved. AnnData
entry points delegate RFU execution to this function and retain attachment
compatibility.

## AnnData Contract

### Reads

- `adata.obsm["airr"]` by default, or another key passed as `airr_key`.
- With `airr_key="TCR_IR"`, a Wells atlas flattened TCR table from
  `adata.uns["TCR_IR"]` or `adata.obsm["TCR_IR"]` (`uns` takes precedence).
- `adata.obs` metadata columns passed to `groupby`, such as `sample`,
  `condition`, or `cell_type`.
- Existing `adata.obs["rfu_label"]` and `adata.obs["rfu_score"]` for summary,
  plotting, and export functions.

### Writes

- `adata.obs["trb_cdr3aa"]`: extracted TRB CDR3 amino-acid sequence.
- `adata.obs["trbv"]`: extracted TRBV/V-gene call.
- `adata.obs["rfu_label"]`: assigned RFU label, missing when unassigned.
- `adata.obs["rfu_score"]`: RFU assignment score, numeric when available.
- `adata.uns["scrfu"]["package_version"]`: installed scRFU package version.
- `adata.uns["scrfu"]["rfu"]`: RFU run provenance by default.

### AIRR Column Aliases

- Cell: `cell_id`, `cell`, `barcode`, `cellid`.
- Chain: `chain`, `locus`.
- CDR3 amino acid: `cdr3aa`, `junction_aa`, `cdr3_aa`, `cdr3`.
- V gene: `v_call`, `v_gene`, `trbv`, `v`.
- Productive: `productive`.

## RFU Calling

### `scrfu.tl.call_rfu`

Purpose: unified RFU calling entrypoint.

Required input: AnnData-like object with AIRR data and `backend="rfu_repo"`.
The RFU checkout is resolved from explicit `rfu_dir`, then `RFU_DIR`; absence of
both is a configuration error.

Main output: pandas DataFrame preserving every extracted row with `input_row_id`,
`unique_sequence_id`, `eligibility_status`, `rfu_id`, `rfu_label`, `rfu_score`,
`pass_thr`, and `rfu_status`; also mutates `adata` by attaching RFU results.

AnnData fields read: `adata.obsm[airr_key]`, or the Wells `TCR_IR` locations
described above, and `adata.obs_names`.

AnnData fields written: `trb_cdr3aa`, `trbv`, `rfu_label`, `rfu_score`, and
`adata.uns["scrfu"]`.

Default mode: `standard`, requiring only official `AssignRFUs()`. Optional
`mode="map_aware"` requires `AssignRFUs_with_map()` and is never selected
implicitly. Exact CDR3 deduplication is enabled by default and can be disabled
with `deduplicate=False`. `chunk_size=None` preserves single-call execution.
Supplying a positive integer enables deterministic serial chunks after
eligibility filtering and deduplication. `resume=True` reuses only fully
validated chunks; `resume=False` reruns them, and `force_recompute=True` takes
precedence over resume.

Common errors: `ValueError` for unknown backend/mode or missing configuration;
`RFUCapabilityError` for an unavailable requested mode; `TypeError` when
`extra_r_args` is a string; downstream extraction errors; `FileNotFoundError`
for an incomplete checkout.

### `scrfu.tl.call_rfu_repo`

Purpose: direct upstream RFU repository backend call.

Required input: AnnData-like object with AIRR data and an external RFU checkout
containing `RFU.R`, `trimerMDSfit_small.Rdata`, and `km5000noMax.Rdata`.

Main output: pandas DataFrame with RFU assignments; attaches results to AnnData.

AnnData fields read: `adata.obsm[airr_key]`, `adata.obs_names`.

AnnData fields written: same as `call_rfu`.

Provenance: stores backend name/mode and resolution source, scRFU version,
resolved directory, detected capabilities, hashes of `RFU.R` and both RFU
`.Rdata` files, wrapper path/hash, timestamp, eligibility rule, deduplication
key, row/query/reconstruction counts, RFU threshold, upstream unique-query miss
count, multiplicity-weighted reconstructed miss count, chain, `airr_key`, and
`out_key`. Chunked runs additionally record the deterministic run ID, chunk
size/count, completed/reused/recomputed/invalidated/failed counts, total elapsed
time, optional maximum resident memory, and run-manifest path.

Common errors: `TypeError` for string `extra_r_args`; `FileNotFoundError` for
missing RFU files or wrapper script; `RuntimeError` if the R backend exits
non-zero or fails to produce output.

### `scrfu.backends.rfu_repo.RFURepoBackend`

Purpose: resolve and validate the external checkout, expose immutable static
capabilities through `backend.capabilities`, enforce explicit backend mode, and
run stable standard or optional map-aware mapping.

Resolution precedence: explicit `rfu_dir`, `RFU_DIR`, actionable configuration
error. `RFURepoPaths.resolution_source` is `explicit_argument` or
`environment_variable`.

Capabilities: immutable booleans for `AssignRFUs`,
`AssignRFUs_with_map`, and `RFUbatch_with_maps`, detected by reading—not
executing—`RFU.R`.

Restartable arguments: `RFURepoBackend.run(..., chunk_size=None, resume=True,
force_recompute=False)`. Chunk sizes must be positive integers. Chunking is
serial and applies to the unique query table, not the original per-cell rows.

### Restartable run layout and cache contract

For chunked execution, `<workdir>/runs/<run_id>/run_manifest.json` indexes
zero-padded directories under `chunks/`. Each chunk directory holds `input.tsv`,
`output.tsv`, `stdout.log`, `stderr.log`, and `manifest.json`. A chunk is reused
only if its manifest parses and is complete; run/chunk identity, offsets, row
count, input and identifier hashes, wrapper/backend/atlas hashes, mode,
threshold, eligibility rule, and deduplication key match; and its output hash,
required columns, row count, unique identifiers, and exact identifier order all
validate. An `output.tsv` by itself is never a cache hit.

The run fingerprint hashes the ordered unique-sequence identifiers and CDR3s,
mode, threshold, deduplication and eligibility rules, RFU and atlas hashes,
wrapper hash/schema, scRFU version, chunk size, and extra R arguments. It omits
timestamps and temporary paths. A chunk ID is
`<run-id-prefix>-<zero-padded-index>-<input-hash-prefix>`. Runtime timestamps do
not affect either identifier.

## Validation and Extraction

### `scrfu.extract.extract_trb_features`

Purpose: extract one TRB CDR3 amino-acid and TRBV value per cell.

Required input: AnnData-like object with a DataFrame-like AIRR table.

Main output: pandas DataFrame with `cell_id`, `cdr3aa`, and `trbv`.

AnnData fields read: `adata.obsm[airr_key]`, or `adata.uns["TCR_IR"]` when
`airr_key="TCR_IR"`, and `adata.obs_names`.

AnnData fields written: none.

Common errors: `KeyError` for missing AIRR key; `ValueError` for missing cell,
chain, CDR3, or V-gene columns.

### `scrfu.extract.extract_wells_tcr_ir_features`

Purpose: adapt a Wells atlas flattened primary VDJ table to the existing scRFU
feature schema without changing RFU result fields.

Required input: `adata.uns["TCR_IR"]` or `adata.obsm["TCR_IR"]` with flattened
VDJ-1 CDR3 and V-call columns, plus `adata.obs_names`.

Main output: pandas DataFrame with `cell_id`, `cdr3aa`, and `trbv`.

AnnData fields written: none.

Common errors: `KeyError` when neither location exists; `ValueError` for missing
flattened fields, ambiguous row alignment, length mismatch, or duplicate cell
identifiers.

### `scrfu.wells.read_wells_receptors_h5ad`

Purpose: read row-aligned Wells `TCR_IR`, observation names, and explicitly
selected metadata without materializing expression or unrelated H5AD elements.

Required input: an AnnData-encoded H5AD containing an encoded `obs` dataframe
and row-aligned `uns["TCR_IR"]` or `obsm["TCR_IR"]` dataframe.

Main output: `WellsReceptorData` containing only `obs`, `tcr_ir`, atlas shape,
container provenance, and source identity.

H5AD fields read: root encoding and `X` shape metadata, selected `obs` fields,
and the selected `TCR_IR` encoded dataframe. `X` values, `raw`, `layers`, `obsp`,
and unrelated `uns` fields are not read.

Common errors: `UnsupportedWellsH5ADLayout` for missing, unsupported, or
misaligned encodings; `ValueError` for invalid selection arguments or metadata
columns.

### `scrfu.wells.prepare_wells_receptor_cache`

Purpose: write checksummed `tcr_ir.tsv.gz`, `obs_metadata.tsv.gz`, and a
fingerprinted `preparation_manifest.json` without copying expression data.

The cache contract and source invalidation semantics are documented in
[wells_receptor_cache.md](wells_receptor_cache.md).

### `scrfu.tl.validate_airr`

Purpose: read-only QC for AIRR-like input before RFU calling.

Required input: AnnData-like object; AIRR data are optional in non-strict mode.

Main output: one-row pandas DataFrame with `airr_key`, `chain`, `n_obs_cells`,
`n_airr_rows`, `n_airr_cells`, `n_overlap_cells`, `barcode_overlap_rate`,
column-presence flags, chain-row counts, productive-chain counts, and `status`.

AnnData fields read: `adata.obsm[airr_key]`, `adata.obs_names` or
`adata.obs.index`.

AnnData fields written: none.

Common statuses: `ok`, `missing_airr_key`, `invalid_airr`,
`missing_required_columns`, `no_chain_rows`, `barcode_mismatch`.

Common errors: only in `strict=True`, raises `ValueError` for a missing AIRR key,
unreadable AIRR table, or missing required columns.

## Attachment

### `scrfu.attach.attach_rfu_results`

Purpose: attach extracted features and RFU assignments back into AnnData.

Required input: AnnData-like object with `.obs` and `.uns`, feature table with
`cell_id` plus CDR3/V columns, and RFU result table with `cell_id`.

Main output: mutates `adata`; returns `None`.

AnnData fields read: `adata.obs.index`.

AnnData fields written: `trb_cdr3aa`, `trbv`, `rfu_label`, `rfu_score`,
`adata.uns["scrfu"]["package_version"]`, and `adata.uns["scrfu"][out_key]`.

Common errors: `ValueError` when feature or RFU tables lack required `cell_id`,
CDR3, or V-gene columns.

## Analysis

Dataset-independent downstream functions are `repertoire_metrics`,
`rfu_metrics`, `rfu_pseudobulk`, `rfu_overlap`,
`rfu_phenotype_coupling`, and `rfu_sequence_matrix`. They require explicit
grouping, sample, and weighting semantics where applicable. See the focused
metric documentation for formulas.

### `scrfu.tl.rfu_summary`

Purpose: summarize RFU assignment rate, unique RFUs, scores, and top RFU.

Required input: AnnData-like object with `adata.obs["rfu_label"]` and
`adata.obs["rfu_score"]`; optional `groupby` column.

Main output: pandas DataFrame, global or grouped.

AnnData fields read: `adata.obs`.

AnnData fields written: none.

Common errors: `ValueError` for missing RFU or grouping columns.

### `scrfu.tl.rfu_metrics`

Purpose: calculate stable descriptive RFU abundance, exact-CDR3 richness,
multiplicity, convergence, entropy, dominance, threshold-pass, and optional
donor/sample prevalence metrics.

Required input: a pandas DataFrame or AnnData-like `.obs`, explicit `groupby`
phenotype columns, and explicit `weighting="cell"` or
`weighting="unique_sequence"`.

Main output: one pandas DataFrame row per phenotype-group/RFU combination.

Input fields read: grouping columns plus configurable cell, CDR3, RFU, threshold,
donor, and sample columns.

Input fields written: none.

Common errors: `TypeError` for unsupported input objects; `ValueError` for
implicit/unknown weighting, empty grouping, or missing columns. Exact metric
definitions are in [rfu_metrics.md](rfu_metrics.md).

### `scrfu.tl.aggregate_rfu`

Purpose: produce a group-by-RFU abundance matrix.

Required input: AnnData-like object with `adata.obs["rfu_label"]` and a grouping
column.

Main output: pandas DataFrame with groups as rows and RFU labels as columns;
values are proportions by default or counts when `normalize=False`.

AnnData fields read: `adata.obs`.

AnnData fields written: none.

Common errors: `ValueError` for missing RFU or grouping columns.

## Plotting

### `scrfu.pl.rfu_bar`

Purpose: bar chart of grouped RFU abundance.

Required input: AnnData-like object accepted by `aggregate_rfu`.

Main output: matplotlib `Axes`.

AnnData fields read: `adata.obs[groupby]`, `adata.obs["rfu_label"]`.

AnnData fields written: none.

Common errors: `ValueError` for missing grouping or RFU columns.

### `scrfu.pl.rfu_heatmap`

Purpose: heatmap of group-by-RFU abundance.

Required input: AnnData-like object accepted by `aggregate_rfu`.

Main output: matplotlib `Axes`.

AnnData fields read: `adata.obs[groupby]`, `adata.obs["rfu_label"]`.

AnnData fields written: none.

Common errors: `ValueError` for missing grouping or RFU columns.

### `scrfu.pl.rfu_score_hist`

Purpose: histogram of non-missing RFU scores.

Required input: AnnData-like object with `adata.obs["rfu_score"]`.

Main output: matplotlib `Axes`.

AnnData fields read: `adata.obs["rfu_score"]`.

AnnData fields written: none.

Common errors: `ValueError` for missing `.obs` or `rfu_score`.

## Export

### `scrfu.io.export_rfu_matrix`

Purpose: write a group-by-RFU matrix to TSV or CSV.

Required input: AnnData-like object accepted by `aggregate_rfu`, grouping column,
and output path.

Main output: pandas DataFrame; writes a delimited file.

AnnData fields read: `adata.obs[groupby]`, `adata.obs["rfu_label"]`.

AnnData fields written: none.

Common errors: `ValueError` for missing grouping or RFU columns; filesystem
errors if the output path cannot be written.
