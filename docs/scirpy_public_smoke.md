# Public Scirpy interoperability smoke

Date: 2026-08-31. This is an interoperability check, not a biological analysis.

Source: Scirpy's openly hosted `wu2020_3k.h5mu` example, a 3,000-cell subset
provided by the Scirpy project. Source URL:
`https://scverse-exampledata.s3.eu-west-1.amazonaws.com/scirpy/wu2020_3k.h5mu`.

| Field | Result |
|---|---:|
| Input size | 17,284,592 bytes |
| Input SHA256 | `a28195a12c9758ac738d3583626694020059d8c33691354ccfa6403110ced566` |
| Cells | 3,000 |
| AIRR chains | 7,544 |
| TRB chains | 3,497 |
| Productive/queryable TRB chains | 2,931 |
| Threshold-qualified TRB chains | 2,246 |
| Native-versus-canonical table rows | 2,931 |
| Native-versus-table mismatches | 0 |
| Fresh native assignment section | 10.52 s |
| Fresh process peak RSS | 607,680 KB |
| Cached assignment section | 1.59 s |
| Cached total command | 3.35 s |
| Cached process peak RSS | 405,924 KB |

The run used Scirpy 0.22.4, scRFU 0.4.0 source on the 0.5 development branch,
standard mode, threshold 0.6, exact-CDR3 deduplication, one worker, and a 5,000
query chunk. It used the untouched official RFU checkout externally. The
written H5MU reloaded successfully, retained all AIRR chains, and passed storage
schema 1.0 validation. The shareable result manifest is named
`official_run/run_manifest.json` in the external smoke output; no downloaded
data or result object is tracked in Git.

