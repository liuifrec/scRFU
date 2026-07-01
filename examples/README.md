# Examples

## Synthetic demo

`synthetic_scirpy_demo.py` builds a tiny synthetic AnnData object with a
scirpy/AIRR-like table, extracts TRB features, attaches synthetic RFU
labels and scores, then runs summary, aggregation, export, and plotting
utilities.

It does not run upstream RFU. The goal is to exercise the AnnData /
scirpy-like data flow, attachment, summary, aggregation, export, and
plotting layers without requiring R, `RFU_DIR`, internet access, or real
10x files.

Run it from the repository root:

```bash
python examples/synthetic_scirpy_demo.py
```

By default it writes these outputs under `.scrfu/demo/`:

- `rfu_matrix.tsv`
- `rfu_bar.png`
- `rfu_heatmap.png`
- `rfu_score_hist.png`

For tests or custom runs, set `SCRFU_DEMO_OUTDIR` to override the output
directory.

## Benchmark h5ad workflow

`benchmark_scirpy_dataset.py` runs manuscript-oriented summaries, matrices,
plots, and a run manifest from a user-provided AnnData/scirpy h5ad file.

```bash
python examples/benchmark_scirpy_dataset.py \
  --input input.h5ad \
  --rfu-dir ~/ext/RFU \
  --groupby sample \
  --cell-type-key cell_type \
  --outdir results/scRFU_benchmark
```

Use `--skip-rfu` only when `adata.obs` already contains `rfu_label` and
`rfu_score`. Use `--dry-run` for input/QC checks without RFU execution.

## GSE190905 preparation scaffold

`gse190905_prepare.py` prepares user-downloaded TCR and metadata CSV/CSV.GZ
files into an h5ad file for later scRFU benchmarking.

```bash
python examples/gse190905_prepare.py \
  --tcr GSE190905_TCR_data.csv.gz \
  --metadata GSE190905_meta_data.csv.gz \
  --out gse190905_scrfu_input.h5ad \
  --report gse190905_prepare_report.json
```

The script infers common column names and supports explicit overrides for TCR
and metadata barcode, chain, CDR3, V-gene, and productive columns.
