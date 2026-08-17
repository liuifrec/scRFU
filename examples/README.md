# Examples

## Dataset-independent receptor workflow

`receptor_table_workflow.py` accepts a canonical receptor cache, canonical TSV,
or H5AD plus an explicit `wells_tcr_ir`/`scirpy_airr` adapter. It writes
canonical receptors, sequence mapping, per-sequence and per-row RFU results,
optional metrics, and a run manifest. `--skip-rfu` supports offline preparation
and output-contract checks without an external RFU checkout.

The Wells workflow remains a case-study wrapper and compatibility output
surface; it now uses the same canonical adapter and table-level RFU engine.

## Offline antigen-evidence workflow

`vdjdb_antigen_evidence.py` uses a user-supplied local VDJdb table with an
explicit release label and optional SHA256. It writes exact long evidence,
ambiguity-preserving summaries, RFU/global coherence, simple grouping
baselines, a size-preserving permutation benchmark, generic figures, and a run
manifest. It does not download a reference or treat a match as proof of antigen
specificity. See [the evidence contract](../docs/antigen_evidence.md).

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
  --rfu-dir /path/to/RFU \
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

## Wells public-atlas workflow

`wells_atlas_workflow.py` accepts either a user-supplied H5AD containing the
Wells `TCR_IR` table or an expression-free prepared receptor cache. Its targeted
HDF5 path reads only selected observation metadata and `TCR_IR`; it does not
materialize `X`, `raw`, layers, or unrelated `uns` fields.

```bash
python examples/wells_atlas_workflow.py \
  --input /path/to/wells_atlas.h5ad \
  --rfu-dir /path/to/RFU \
  --outdir results/wells_scrfu \
  --chunk-size 20000 \
  --primary-chain
```

Prepare a reusable cache first when the same receptor data will be analyzed
repeatedly:

```bash
scrfu prepare-wells /path/to/wells_atlas.h5ad \
  --output-dir /path/to/wells_receptor_cache \
  --obs-column donor_id \
  --obs-column library_id \
  --obs-column cell_type

python examples/wells_atlas_workflow.py \
  --input /path/to/wells_receptor_cache \
  --source-h5ad /path/to/wells_atlas.h5ad \
  --rfu-dir /path/to/RFU \
  --outdir results/wells_scrfu \
  --chunk-size 20000
```

The optional `--source-h5ad` verifies that the cache fingerprint still matches
the source. See [the cache design](../docs/wells_receptor_cache.md).

`--rfu-dir` may be omitted when `RFU_DIR` is set. Use `--max-cells 1000` for a
small integration smoke test, `--no-resume` to ignore completed chunks, or
`--force-recompute` to force all chunks to be replaced. `--all-productive-chains`
retains productive TRB rows from both VDJ slots; its per-cell output can contain
more than one row for a cell. `--write-annotated` is currently limited to the
primary-chain workflow.

The workflow writes `extracted_trb.tsv.gz`, `adapter_qc.json`,
`unique_sequence_map.tsv.gz` (one mapping row per extracted input row),
sequence- and cell-level RFU result tables, and a top-level `run_manifest.json`.
Backend chunk artifacts are retained below `<outdir>/backend/runs/<run_id>/`.
`--write-annotated` requires original H5AD input and explicitly opts into
expression loading; a receptor cache cannot produce an annotated H5AD.

A 1,000-cell smoke test and a complete restartable run differ only by the
optional cell limit:

```bash
RFU_DIR=/path/to/RFU python examples/wells_atlas_workflow.py \
  --input /path/to/wells_atlas.h5ad \
  --outdir results/wells_scrfu_smoke \
  --chunk-size 20000 \
  --max-cells 1000 \
  --primary-chain

RFU_DIR=/path/to/RFU python examples/wells_atlas_workflow.py \
  --input /path/to/wells_atlas.h5ad \
  --outdir results/wells_scrfu \
  --chunk-size 20000 \
  --primary-chain
```
