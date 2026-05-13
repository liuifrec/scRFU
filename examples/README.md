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
