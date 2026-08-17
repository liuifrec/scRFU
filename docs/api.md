# API Reference

The intended stable public API is documented in
[api_contract.md](api_contract.md).

Current public modules:

- `scrfu.pp`: canonical receptor schema and validation.
- `scrfu.adapters`: named receptor adapters and optional modality routing.
- `scrfu.tl`: RFU calling, validation, downstream metrics, antigen evidence,
  longitudinal methods, transfer diagnostics, robustness, and comparators.
- `scrfu.backends.rfu_repo`: backend path resolution, immutable capability
  detection, explicit mode enforcement, and stable RFU mapping.
- `scrfu.pl`: matplotlib plotting helpers.
- `scrfu.io`: h5ad helpers and RFU matrix export.

scRFU calls a user-provided upstream RFU repository through
`backend="rfu_repo"` and does not vendor upstream RFU code or data.

The exact public upstream behavior, optional map-aware capability boundary, and
Wells atlas adapter are documented in
[upstream_rfu_semantics.md](upstream_rfu_semantics.md).
Memory-efficient Wells extraction and cache validation are documented in
[wells_receptor_cache.md](wells_receptor_cache.md). Stable descriptive RFU
metrics and their explicit grouping and weighting definitions are documented in
[rfu_metrics.md](rfu_metrics.md).

The complete stable/experimental/compatibility classification is in
[api_freeze_v1.md](api_freeze_v1.md). Longitudinal and frozen-transfer APIs are
experimental pending independent real-cohort validation.
