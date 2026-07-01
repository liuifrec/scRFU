# scRFU

scRFU is an AnnData/scirpy-compatible Python framework for applying established
RFU-based TCR functional annotation inside single-cell immune profiling
workflows.

## What scRFU Does

- Extracts TRB CDR3 amino-acid and TRBV features from AnnData/scirpy-style AIRR
  tables.
- Calls a user-provided upstream RFU repository backend.
- Attaches RFU labels, RFU scores, extracted TRB features, and provenance into
  AnnData.
- Provides RFU summaries, group-level RFU matrices, matplotlib plots, and TSV
  export for downstream analysis.
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

## External RFU Dependency

RFU assignment requires a separate checkout of the upstream RFU repository. Pass
that location with `rfu_dir` or the example scripts' `--rfu-dir` option. scRFU
records hashes of upstream RFU files when RFU is run, but the files remain
outside this repository.

## Quickstart

```python
import scrfu

qc = scrfu.tl.validate_airr(adata, airr_key="airr", chain="TRB")

scrfu.tl.call_rfu(
    adata,
    backend="rfu_repo",
    rfu_dir="~/ext/RFU",
)

summary = scrfu.tl.rfu_summary(adata, groupby="sample")
matrix = scrfu.tl.aggregate_rfu(adata, groupby="sample")

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

## Analysis Utilities

The public analysis layer includes `scrfu.tl.rfu_summary`,
`scrfu.tl.aggregate_rfu`, `scrfu.pl.rfu_bar`, `scrfu.pl.rfu_heatmap`,
`scrfu.pl.rfu_score_hist`, and `scrfu.io.export_rfu_matrix`.

## Reproducible Examples

- [examples/synthetic_scirpy_demo.py](examples/synthetic_scirpy_demo.py):
  offline synthetic AnnData/scirpy-style demo.
- [examples/real_scirpy_workflow.py](examples/real_scirpy_workflow.py):
  template for user-provided real h5ad input.
- [examples/radiation_pbmc_workflow.py](examples/radiation_pbmc_workflow.py):
  radiation-associated PBMC workflow template without bundled data.

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
