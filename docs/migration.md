# RFU Result Migration Note

## Headerless upstream input

Older scRFU wrappers could write a `CDR3` column header into the temporary file
passed to upstream `EncodeRepertoire()`. Upstream reads that file with
`header=FALSE`, and the text `CDR3` starts with `C`, so the header could be
encoded as a synthetic repertoire row. The current wrapper writes the upstream
file without column names. RFU results made with an older scRFU wrapper should
therefore be regenerated rather than mixed with current output.

The current standard path also deduplicates eligible queries by exact `cdr3aa`
before calling RFU. Integration tests against an official-compatible checkout
showed identical nearest RFU IDs and labels, threshold-pass flags, and scores to
absolute tolerance `1e-12`, while reconstruction restores the original
cell-level multiplicity and row order. Rows whose CDR3 does not start with `C`
are retained with explicit eligibility and RFU status instead of silently
disappearing.

These are deliberate behavior and public-API changes. This repository has no
formal versioning policy beyond its current pre-alpha `0.x` status. A minor
version increment from `0.1.x` to `0.2.0` is recommended for the release that
includes them; this documentation change does not itself change the package
version.
