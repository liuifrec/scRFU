# Scirpy/AnnData/MuData workflow

Start with the [synthetic native scverse tutorial](../examples/tutorial_scverse_native.py):

```bash
python -m pip install ".[scirpy,mudata,plotting]"
python examples/tutorial_scverse_native.py --outdir /tmp/scrfu-scverse-tutorial
```

This executable example constructs synthetic AIRR and MuData objects. Its mock
backend and `mock_RFU` labels test chain alignment, storage and cell summaries;
they are not scientific assignments. For actual assignment, configure the
external reference as described in [first real RFU run](first_real_rfu.md).
Study notebooks are maintained separately from this software distribution.

```python
import scanpy as sc
import scirpy as ir
import scrfu

# `mdata` contains an AIRR modality produced by Scirpy.
# Chain indexing is optional for assignment; use it for Scirpy QC if needed.
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
columns belong to the AIRR modality. Transfer them explicitly to the GEX
observation index before using its existing RNA embedding:

```python
airr = mdata.mod["airr"]
gex = mdata.mod["gex"]
columns = list(airr.obs.filter(like="scrfu_"))
for column in columns:
    gex.obs[column] = airr.obs[column].reindex(gex.obs_names)
# Choose the embedding that actually exists in the GEX object.
sc.pl.embedding(gex, basis="umap_orig", color="scrfu_assignment_status")
mdata.pull_obs(mods=["airr"], columns=columns)
```

This preserves RNA-only rows instead of using receptor-bearing cells as the RNA
denominator. The API distinguishes absent, noneligible, threshold-miss and
ambiguous cells, and supports an explicit cell-summary policy before RFU pseudobulk. A nearest
RFU label alone does not establish qualification; the score is not calibrated
confidence. Request metadata columns explicitly after annotation; do not rely on implicit
MuData metadata propagation.

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
