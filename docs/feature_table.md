# Feature Comparison

This table compares the original upstream RFU workflow with scRFU. scRFU
depends on the upstream RFU method for RFU assignment and does not replace it.

| Feature | Upstream RFU R workflow | scRFU |
| --- | --- | --- |
| RFU algorithm | Original RFU algorithm and atlas | Uses upstream RFU for assignment |
| Requires external RFU files | Yes | Yes, via `rfu_dir` |
| AnnData integration | Not the primary interface | Yes |
| scirpy/AIRR extraction | Not the primary interface | Yes, AnnData/scirpy-style AIRR input |
| Per-cell RFU annotation | RFU output can be mapped externally | Yes, stored in `adata.obs` |
| Group-level RFU matrix | Requires downstream handling | Yes |
| Plotting | RFU-specific workflow outputs | Lightweight matplotlib helpers |
| Export | RFU workflow outputs | Group-by-RFU matrix export |
| Provenance tracking | External workflow dependent | Stored in `adata.uns["scrfu"]` |
| Synthetic demo | Not the primary package focus | Yes, CI-safe synthetic demo |
| Real-data workflow template | User-managed workflow | Yes, user-provided input paths |
