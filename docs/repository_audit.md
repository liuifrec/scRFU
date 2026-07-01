# Repository Audit

## Current Implemented Features

- External upstream RFU repository backend integration through
  `scrfu.tl.call_rfu` and `scrfu.tl.call_rfu_repo`.
- AnnData/scirpy-like AIRR extraction for TRB CDR3 amino-acid and TRBV features.
- Read-only AIRR input QC through `scrfu.tl.validate_airr`.
- RFU assignment attachment into `adata.obs` with provenance under
  `adata.uns["scrfu"]`.
- RFU assignment summaries, group-level aggregation, matplotlib plots, and RFU
  matrix export.
- Offline synthetic scirpy-style demo.
- Real-data workflow templates and benchmark scaffolds for user-provided h5ad or
  CSV/CSV.GZ inputs.
- Planning documentation for manuscript, roadmap, ALFU future design, benchmark
  outputs, journal positioning, and figure generation.

## Public API Surface

- RFU calling: `scrfu.tl.call_rfu`, `scrfu.tl.call_rfu_repo`.
- Validation and extraction: `scrfu.tl.validate_airr`,
  `scrfu.extract.extract_trb_features`.
- Attachment: `scrfu.attach.attach_rfu_results`.
- Analysis: `scrfu.tl.rfu_summary`, `scrfu.tl.aggregate_rfu`.
- Plotting: `scrfu.pl.rfu_bar`, `scrfu.pl.rfu_heatmap`,
  `scrfu.pl.rfu_score_hist`.
- Export: `scrfu.io.export_rfu_matrix`.
- CLI: `scrfu call-rfu`.
- Example scripts are intended as reproducible workflow scaffolds, not as a
  stable library API.

## Current Tests and Coverage Areas

- AIRR extraction aliases, missing-key behavior, and TRB filtering.
- Read-only AIRR validation statuses for valid input, missing keys, missing
  required columns, absent TRB rows, barcode mismatch, and strict mode.
- RFU attachment behavior, alternate column names, empty features, and
  provenance storage.
- RFU backend path validation and provenance hash/timestamp fields.
- RFU summary, aggregation, plotting, and matrix export.
- CLI argument parsing and dispatch.
- Synthetic demo execution and output creation.
- Workflow script help/error paths.
- Benchmark script dry-run and skip-RFU workflows using tiny synthetic h5ad
  files.
- GSE190905 preparation using tiny offline gzip CSV inputs.

## Known Limitations

- scRFU depends on a user-provided upstream RFU checkout and does not validate
  upstream licensing or installation beyond required file presence.
- RFU scientific behavior is inherited from the upstream RFU method; scRFU does
  not add a new RFU model or atlas.
- AnnData AIRR handling is intentionally conservative and optimized for
  DataFrame-like scirpy/AIRR tables used by the current workflows.
- The GSE190905 preparation scaffold stores one representative AIRR row per cell
  for current scRFU input compatibility; raw multi-contig details are not
  preserved in `obsm`.
- Real-data biological demonstrations are not completed in-repository.
- No runtime/resource benchmark has been generated on a real dataset yet.
- No comparison to clonotype, diversity, V/J usage, or other conventional
  repertoire metrics is implemented yet.
- Package metadata still indicates pre-alpha status and placeholder repository
  URLs.

## Missing Items for Independent Manuscript Submission

- At least two completed real-data demonstrations with documented inputs,
  manifests, and generated figures.
- A conventional repertoire-metric comparison layer or companion analysis
  scripts.
- Runtime and memory summaries across synthetic and real datasets.
- Versioned release artifacts and stable package metadata.
- A clear licensing statement for interaction with the upstream RFU dependency.
- Reproducible manuscript figure scripts that produce final panels from saved
  benchmark outputs.
- Review of edge cases for full scirpy awkward-array AIRR inputs.
- Broader documentation examples for common metadata layouts and failure modes.

## Recommended Pre-Submission Checklist

- Run the synthetic demo and keep outputs reproducible.
- Run `benchmark_scirpy_dataset.py` on a public COVID-19 PBMC single-cell TCR
  dataset.
- Download and inspect GSE190905 processed TCR and metadata files, then run
  `gse190905_prepare.py` and the benchmark workflow.
- Compare RFU matrices with clonotype abundance, diversity, V gene usage, and
  cell-type composition.
- Record runtime, memory, RFU assignment rate, and input QC tables.
- Generate all manuscript figures from scripts and archived result files.
- Confirm upstream RFU citation and licensing language.
- Freeze a v0.x release, update package URLs/authors, and tag the exact code
  used for manuscript figures.
