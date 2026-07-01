# API Contract

This document describes the intended stable public API for manuscript reviewers
and future users. scRFU is currently pre-alpha, but these functions are the
planned compatibility surface.

## AnnData Contract

### Reads

- `adata.obsm["airr"]` by default, or another key passed as `airr_key`.
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

Required input: AnnData-like object with AIRR data, `backend="rfu_repo"`, and
`rfu_dir` pointing to a user-provided upstream RFU checkout.

Main output: pandas DataFrame with at least `cell_id`, `rfu_label`, and
`rfu_score`; also mutates `adata` by attaching RFU results.

AnnData fields read: `adata.obsm[airr_key]`, `adata.obs_names`.

AnnData fields written: `trb_cdr3aa`, `trbv`, `rfu_label`, `rfu_score`, and
`adata.uns["scrfu"]`.

Common errors: `ValueError` for unknown backend or missing `rfu_dir`; `TypeError`
when `extra_r_args` is a string; downstream `KeyError` or `ValueError` from AIRR
extraction; `FileNotFoundError` for incomplete upstream RFU checkout.

### `scrfu.tl.call_rfu_repo`

Purpose: direct upstream RFU repository backend call.

Required input: AnnData-like object with AIRR data and an external RFU checkout
containing `RFU.R`, `trimerMDSfit_small.Rdata`, and `km5000noMax.Rdata`.

Main output: pandas DataFrame with RFU assignments; attaches results to AnnData.

AnnData fields read: `adata.obsm[airr_key]`, `adata.obs_names`.

AnnData fields written: same as `call_rfu`.

Provenance: stores backend name, scRFU version, upstream RFU directory, hashes of
`RFU.R` and RFU `.Rdata` files, wrapper R script path/hash, timestamp, chain,
`airr_key`, and `out_key`.

Common errors: `TypeError` for string `extra_r_args`; `FileNotFoundError` for
missing RFU files or wrapper script; `RuntimeError` if the R backend exits
non-zero or fails to produce output.

## Validation and Extraction

### `scrfu.extract.extract_trb_features`

Purpose: extract one TRB CDR3 amino-acid and TRBV value per cell.

Required input: AnnData-like object with a DataFrame-like AIRR table.

Main output: pandas DataFrame with `cell_id`, `cdr3aa`, and `trbv`.

AnnData fields read: `adata.obsm[airr_key]`, `adata.obs_names`.

AnnData fields written: none.

Common errors: `KeyError` for missing AIRR key; `ValueError` for missing cell,
chain, CDR3, or V-gene columns.

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

### `scrfu.tl.rfu_summary`

Purpose: summarize RFU assignment rate, unique RFUs, scores, and top RFU.

Required input: AnnData-like object with `adata.obs["rfu_label"]` and
`adata.obs["rfu_score"]`; optional `groupby` column.

Main output: pandas DataFrame, global or grouped.

AnnData fields read: `adata.obs`.

AnnData fields written: none.

Common errors: `ValueError` for missing RFU or grouping columns.

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
