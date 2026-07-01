# Manuscript Figure Plan

This plan avoids claiming biological findings before real-data analyses are
completed.

## Figure 1: Software Architecture and Workflow

- Required input: conceptual schematic plus documented API contract.
- Script to generate: future manuscript figure script or vector graphics source.
- Expected output files: architecture figure panel, API/provenance summary.
- Current status: needs figure artwork; software components are implemented.

Content: AnnData/scirpy AIRR input, TRB feature extraction, external upstream RFU
execution, RFU attachment into `adata.obs`, provenance in `adata.uns`, grouped
summary, visualization, export, and downstream scverse-compatible analysis.

## Figure 2: Synthetic/scirpy Demo Validation

- Required input: synthetic AnnData object from
  `examples/synthetic_scirpy_demo.py`.
- Script to generate: `python examples/synthetic_scirpy_demo.py`.
- Expected output files: `rfu_matrix.tsv`, `rfu_bar.png`, `rfu_heatmap.png`,
  `rfu_score_hist.png`.
- Current status: implemented with synthetic RFU labels; does not run upstream
  RFU.

Content: feature extraction, RFU-style annotation attachment, summary,
aggregation, plotting, and export on a CI-safe demo.

## Figure 3: COVID-19 PBMC Benchmark

- Required input: user-provided public COVID-19 PBMC AnnData/scirpy h5ad with
  TCR AIRR data and metadata such as disease, severity, sample, and cell type.
- Script to generate:
  `python examples/benchmark_scirpy_dataset.py --input input.h5ad --rfu-dir ~/ext/RFU --groupby sample --cell-type-key cell_type --outdir results/scRFU_benchmark`.
- Expected output files: `validate_airr.tsv`, `rfu_summary_by_group.tsv`,
  `rfu_matrix_by_group.tsv`, `rfu_bar_by_group.png`,
  `rfu_heatmap_by_group.png`, `rfu_score_hist.png`, optional cell-type outputs,
  and `run_manifest.json`.
- Current status: workflow scaffold implemented; needs real data and manuscript
  panel assembly.

Content: RFU assignment rate, RFU composition by disease/severity/sample when
metadata permit, RFU score distribution, and optional cell-type association.

## Figure 4: Radiotherapy-Associated Public Dataset GSE190905

- Required input: user-downloaded processed GSE190905 TCR and metadata files,
  then prepared h5ad.
- Script to generate:
  `python examples/gse190905_prepare.py --tcr GSE190905_TCR_data.csv.gz --metadata GSE190905_meta_data.csv.gz --out gse190905_scrfu_input.h5ad --report gse190905_prepare_report.json`,
  followed by `benchmark_scirpy_dataset.py` or
  `examples/radiation_pbmc_workflow.py`.
- Expected output files: preparation report JSON, prepared h5ad, AIRR QC,
  grouped RFU summaries, RFU matrices, heatmaps, score histogram, and manifest.
- Current status: preparation scaffold implemented; needs real downloaded files,
  metadata inspection, RFU run, and analysis validation.

Content: input QC, pre/post or treatment-group RFU composition if supported by
metadata, sample-level RFU heatmap, and assignment summary.

## Supplementary Figure: API, Provenance, and Output Schema

- Required input: `docs/api_contract.md`, `docs/benchmark_outputs.md`, and saved
  run manifests.
- Script to generate: future manuscript supplement script.
- Expected output files: API/provenance diagram, output schema table, example
  manifest excerpt.
- Current status: documentation implemented; figure/table generation remains
  future work.
