# API Reference

The intended stable public API is documented in
[api_contract.md](api_contract.md).

Current public modules:

- `scrfu.tl`: RFU calling, AIRR validation, summary, and aggregation.
- `scrfu.extract`: AIRR/TRB feature extraction.
- `scrfu.attach`: RFU result attachment into AnnData.
- `scrfu.pl`: matplotlib plotting helpers.
- `scrfu.io`: h5ad helpers and RFU matrix export.

scRFU calls a user-provided upstream RFU repository through
`backend="rfu_repo"` and does not vendor upstream RFU code or data.
