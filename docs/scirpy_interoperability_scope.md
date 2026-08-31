# Scirpy interoperability and scope

Scirpy is the general immune-repertoire ecosystem layer. It imports AIRR data,
stores every receptor chain in `adata.obsm["airr"]`, indexes receptor-chain
models, defines clonotypes and distances, and provides broad repertoire
analysis and visualization. scRFU does not replace those capabilities.

scRFU adds a frozen-reference, exact-CDR3-based TRB representation and
RFU-specific analyses. `scrfu.tl.assign_rfu` accepts the current Scirpy Awkward
AIRR representation directly and writes one aligned result record per AIRR
chain. It delegates every scientific assignment to the canonical table-level
RFU implementation.

| API area | Classification | Recommended use |
|---|---|---|
| RFU assignment, coverage and frozen-reference provenance | Core RFU method | scRFU |
| RFU pseudobulk, convergence, overlap and phenotype coupling | RFU-specific downstream | scRFU |
| AIRR Awkward extraction, chain alignment and MuData routing | Interoperability helper | scRFU with Scirpy objects |
| Repertoire diversity and generic summaries | Generic repertoire convenience | Scirpy for broad workflows; retained in scRFU for stable compatibility |
| Historical one-row-per-cell AnnData wrapper | Compatibility API | New work should prefer `assign_rfu` |

`scirpy.pp.index_chains` is not required for chain-level RFU assignment. If
present, its `obsm["chain_indices"]` output is consulted only when the user
explicitly requests the `primary_vdj` cell-summary policy. Running or not
running `index_chains` cannot change chain-level RFU results.

The external official RFU implementation remains required for real
assignments. Scirpy, AnnData and MuData do not redistribute those assets.

