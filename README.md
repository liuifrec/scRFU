# scRFU

scRFU is a receptor-table and AnnData/scirpy-compatible Python framework for
applying established RFU-based TCR functional annotation inside single-cell
immune profiling workflows. Its core is dataset-independent; Wells is a named
adapter, benchmark, and reproducible case study rather than the data model.

## What scRFU Does

- Extracts TRB CDR3 amino-acid and TRBV features from AnnData/scirpy-style AIRR
  tables.
- Prepares canonical receptor rows through Wells, scirpy/AIRR, generic
  DataFrame, and Cell Ranger VDJ adapters while keeping metadata separate.
- Selectively reads receptor and observation data without loading expression
  matrices and writes portable checksummed receptor caches.
- Calls a user-provided upstream RFU repository backend.
- Attaches RFU labels, RFU scores, extracted TRB features, and provenance into
  AnnData.
- Provides RFU summaries, group-level RFU matrices, matplotlib plots, and TSV
  export for downstream analysis.
- Provides conventional repertoire metrics, sample-level RFU pseudobulk,
  overlap, phenotype coupling, explicit convergence definitions, and generic
  manuscript-oriented plots.
- Provides CI-safe synthetic examples plus real-data workflow scaffolds for
  user-provided datasets.

## What scRFU Does Not Do

- scRFU does not invent RFU or replace the upstream RFU method.
- scRFU does not vendor upstream RFU code, atlas files, or `.Rdata` objects.
- scRFU does not download real datasets.
- B-cell/ALFU support is design-only future work.
- Radiation biology claims require completed real-data analyses and are not
  implied by the workflow templates.

## Installation

Development install:

```bash
pip install -e ".[dev]"
```

Optional scirpy extras:

```bash
pip install -e ".[dev,scirpy]"
```

Core imports do not require matplotlib. Install plotting support with
`pip install -e ".[plotting]"`; `.[all]` installs plotting and scirpy extras.

## External RFU Dependency

RFU assignment requires a separate checkout of the upstream RFU repository.
Pass that location with `rfu_dir`/`--rfu-dir`, or set `RFU_DIR`. An explicit
argument takes precedence over the environment variable. Standard mode uses the
official public `AssignRFUs()` implementation and remains the default. scRFU
records hashes and detected optional capabilities, but the files remain outside
this repository.

## Quickstart

```python
import scrfu

prepared = scrfu.adapters.prepare_receptors(
    airr_dataframe,
    adapter="generic_airr_dataframe",
    chain="TRB",
)
table_run = scrfu.tl.call_rfu_table(
    prepared.receptors,
    rfu_dir="/path/to/RFU",
    chunk_size=20_000,
)

qc = scrfu.tl.validate_airr(adata, airr_key="airr", chain="TRB")

scrfu.tl.call_rfu(
    adata,
    backend="rfu_repo",
    rfu_dir="/path/to/RFU",
    mode="standard",
    chunk_size=20_000,
    workdir="results/rfu_work",
)

summary = scrfu.tl.rfu_summary(adata, groupby="sample")
matrix = scrfu.tl.aggregate_rfu(adata, groupby="sample")

pseudobulk = scrfu.tl.rfu_pseudobulk(
    per_row_results,
    sample_key="donor_id",
    weighting="cell",
    normalize="proportion",
)
overlap = scrfu.tl.rfu_overlap(pseudobulk, metric="jaccard")

scrfu.pl.rfu_bar(adata, groupby="sample")
scrfu.pl.rfu_heatmap(adata, groupby="sample")
scrfu.pl.rfu_score_hist(adata)
```

## Input Schema

scRFU reads AIRR-like data from `adata.obsm["airr"]`. Supported aliases include
cell barcode (`cell_id`, `cell`, `barcode`, `cellid`), chain (`chain`, `locus`),
CDR3 amino-acid sequence (`cdr3aa`, `junction_aa`, `cdr3_aa`, `cdr3`), V gene
(`v_call`, `v_gene`, `trbv`, `v`), and optional `productive`.

RFU results are written to:

- `adata.obs["trb_cdr3aa"]`
- `adata.obs["trbv"]`
- `adata.obs["rfu_label"]`
- `adata.obs["rfu_score"]`
- `adata.uns["scrfu"]`

See [docs/api_contract.md](docs/api_contract.md) for the stable API and AnnData
contract.

The canonical schema, adapters, and expression-free cache are documented in
[docs/receptor_schema.md](docs/receptor_schema.md),
[docs/adapters.md](docs/adapters.md), and
[docs/receptor_cache.md](docs/receptor_cache.md). Expression matrices are not
required after receptor preparation.

## Analysis Utilities

The public analysis layer includes legacy summaries plus
`repertoire_metrics`, `rfu_metrics`, `rfu_pseudobulk`, `rfu_overlap`,
`rfu_phenotype_coupling`, and generic `scrfu.pl` functions. See
[docs/v0.2_api_review.md](docs/v0.2_api_review.md) for the compatibility
boundary and the focused metric documentation for formulas.

## Reproducible Examples

- [examples/synthetic_scirpy_demo.py](examples/synthetic_scirpy_demo.py):
  offline synthetic AnnData/scirpy-style demo.
- [examples/real_scirpy_workflow.py](examples/real_scirpy_workflow.py):
  template for user-provided real h5ad input.
- [examples/radiation_pbmc_workflow.py](examples/radiation_pbmc_workflow.py):
  radiation-associated PBMC workflow template without bundled data.
- [examples/wells_atlas_workflow.py](examples/wells_atlas_workflow.py):
  restartable Wells public-atlas workflow for a user-supplied H5AD or compact
  expression-free receptor cache.
- [examples/receptor_table_workflow.py](examples/receptor_table_workflow.py):
  dataset-independent cache, canonical TSV, or named-H5AD-adapter workflow.
- [examples/rfu_downstream_analysis.py](examples/rfu_downstream_analysis.py):
  sample-level downstream metrics and generic figures from RFU per-row results.
- [examples/benchmark_summary.py](examples/benchmark_summary.py): combines
  user-labelled manifests into a manuscript-ready benchmark table.

For large inputs, chunking occurs after exact-CDR3 deduplication. Completed
chunks are reused only after their manifest, scientific inputs, hashes, output
schema, row count, and identifiers all validate. `chunk_size=None` retains the
existing one-call behavior; `resume=False` reruns every chunk, and
`force_recompute=True` takes precedence and forces every chunk to run again.
Chunk execution is currently serial.

## Manuscript and Benchmark Workflows

- [examples/benchmark_scirpy_dataset.py](examples/benchmark_scirpy_dataset.py):
  benchmark scaffold for user-provided AnnData/scirpy h5ad files.
- [examples/gse190905_prepare.py](examples/gse190905_prepare.py):
  preparation scaffold for user-downloaded GSE190905-like TCR and metadata
  files.
- [docs/benchmark_outputs.md](docs/benchmark_outputs.md): output schema.
- [docs/repository_audit.md](docs/repository_audit.md): manuscript-readiness
  audit.
- [docs/manuscript_figures.md](docs/manuscript_figures.md): figure-generation
  plan.

## Citation and Attribution

Please cite and attribute the original upstream RFU method when using RFU
assignments. scRFU is an integration framework around user-provided upstream RFU
code/data and should be cited separately once a release or manuscript is
available.
