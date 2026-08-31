# Scirpy/AnnData/MuData workflow

The runnable source is
[`examples/tutorial_scverse_native.py`](../examples/tutorial_scverse_native.py).
Its default results are deterministic synthetic assignments for software
testing, not official RFU results.

```python
import scanpy as sc
import scirpy as ir
import scrfu

# `mdata` contains an AIRR modality produced by Scirpy.
ir.pp.index_chains(mdata)
scrfu.tl.assign_rfu(
    mdata,
    airr_mod="airr",
    airr_key="airr",
    key_added="scrfu",
    cell_summary=True,
    summary_policy="ambiguity_aware",
    rfu_dir="/user/supplied/upstream/RFU",
)
```

Every source chain is still present in
`mdata.mod["airr"].obsm["airr"]`. The aligned result is in
`mdata.mod["airr"].obsm["scrfu"]`. TRA, BCR and nonproductive chains have
explicit non-eligible records; multiple TRBs are not collapsed.

Cell summaries are opt-in because a cell can contain conflicting TRB
assignments. Once requested, the `scrfu_rfu_label` and `scrfu_rfu_score`
columns can be grouped or plotted normally:

```python
airr = mdata.mod["airr"]
sc.pl.umap(airr, color="scrfu_rfu_label")
airr.obs.groupby(["cell_type", "scrfu_rfu_label"], observed=True).size()
```

RFU pseudobulk, phenotype coupling, longitudinal utilities and reference
coverage operate on explicit receptor/result tables; the tutorial demonstrates
the conversion without accessing expression values. Save/reload with ordinary
`write_h5ad`/`read_h5ad` or `write_h5mu`/`read_h5mu`, then call
`scrfu.tl.validate_scrfu_schema`.

Install `.[scirpy,mudata,plotting]` for this workflow. Real RFU assignment
requires the separately obtained official RFU checkout and its R runtime. The
package does not redistribute those assets. The targeted receptor-only H5AD
reader remains the supported path for very large atlases; native backed
mutation is not currently claimed.

