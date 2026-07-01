# Benchmark Output Schema

This document describes expected outputs from manuscript-oriented scRFU
workflow scripts. The scripts do not download data and should be run on
user-provided inputs.

## `benchmark_scirpy_dataset.py`

| Filename | Purpose | Expected columns or contents | Manuscript use |
| --- | --- | --- | --- |
| `validate_airr.tsv` | Input AIRR QC table | `airr_key`, `chain`, `n_obs_cells`, `n_airr_rows`, `n_airr_cells`, `n_overlap_cells`, `barcode_overlap_rate`, column flags, `n_chain_rows`, `n_productive_chain_rows`, `status` | Input QC and Supplementary methods |
| `rfu_summary_global.tsv` | Dataset-level assignment summary | `n_cells`, `n_assigned`, `assignment_rate`, `n_unique_rfu`, score summaries, `top_rfu`, `top_rfu_count` | Assignment-rate table |
| `rfu_summary_by_group.tsv` | Grouped assignment summary | Group column plus global summary columns | Figure 3 or Figure 4 assignment panels |
| `rfu_matrix_by_group.tsv` | Group-by-RFU abundance matrix | First column is group index; RFU labels as remaining columns | Heatmap input |
| `rfu_bar_by_group.png` | Grouped RFU abundance plot | Image file | Composition panel |
| `rfu_heatmap_by_group.png` | Group-by-RFU heatmap | Image file | Main heatmap panel |
| `rfu_score_hist.png` | RFU score distribution | Image file | Score/QC panel |
| `rfu_summary_by_cell_type.tsv` | Optional cell-type assignment summary | Cell-type column plus summary columns | Optional cell-state association |
| `rfu_matrix_by_cell_type.tsv` | Optional cell-type RFU matrix | Cell types by RFUs | Optional heatmap |
| `rfu_heatmap_by_cell_type.png` | Optional cell-type heatmap | Image file | Optional cell-state panel |
| `run_manifest.json` | Reproducibility metadata | scRFU version, input path, output directory, RFU directory, RFU run flag, grouping keys, cell count, output files, timestamp, external RFU note | Reproducibility and supplement |
| `annotated.h5ad` | Optional annotated AnnData | AnnData with RFU columns and provenance | Reanalysis artifact |

## `gse190905_prepare.py`

| Filename | Purpose | Expected columns or contents | Manuscript use |
| --- | --- | --- | --- |
| User-selected `.h5ad` output | Prepared AnnData input for scRFU | `adata.obs` metadata; row-aligned `adata.obsm["airr"]` with `cell_id`, `chain`, `cdr3aa`, `v_call`, optional `productive` | Figure 4 input |
| User-selected report JSON | Preparation QC report | `n_metadata_rows`, `n_tcr_rows`, `n_cells_with_tcr`, `n_overlap_cells`, `barcode_overlap_rate`, `inferred_columns`, `n_productive_trb_rows`, `output_path` | Input QC table and methods |

The prepared AnnData stores one representative AIRR row per cell for current
scRFU compatibility. The report records raw TCR row counts before this
row-aligned representation is created.

## `radiation_pbmc_workflow.py`

| Filename | Purpose | Expected columns or contents | Manuscript use |
| --- | --- | --- | --- |
| `rfu_summary_by_sample.tsv` | Sample-level RFU assignment summary | `sample` plus assignment summary columns | Figure 4 sample summary |
| `rfu_matrix_by_sample.tsv` | Sample-by-RFU matrix | Sample index plus RFU columns | Sample-level heatmap |
| `rfu_summary_by_condition.tsv` | Condition-level RFU assignment summary | `condition` plus assignment summary columns | Treatment/time/dose summary |
| `rfu_matrix_by_condition.tsv` | Condition-by-RFU matrix | Condition index plus RFU columns | Treatment/time/dose heatmap |
| `rfu_bar_by_sample.png` | Sample-level RFU bar plot | Image file | Exploratory figure |
| `rfu_heatmap_by_sample.png` | Sample-level RFU heatmap | Image file | Figure 4 heatmap candidate |
| `rfu_bar_by_condition.png` | Condition-level RFU bar plot | Image file | Exploratory condition panel |
| `rfu_heatmap_by_condition.png` | Condition-level RFU heatmap | Image file | Figure 4 condition panel |
| `rfu_score_hist.png` | RFU score distribution | Image file | QC panel |
| `annotated.h5ad` | Annotated AnnData | RFU columns and provenance | Reanalysis artifact |
| `rfu_work/` | RFU backend scratch files and logs | Input/output TSVs and R logs when produced | Debugging and reproducibility |

The radiation workflow is currently a template. It does not by itself establish
radiation-associated biological findings.
