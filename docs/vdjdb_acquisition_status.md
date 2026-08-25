# VDJdb acquisition status

Updated on 2026-08-25. The explicitly pinned public release was downloaded to
an external runtime directory, validated, and analyzed. It remains absent from
the repository, wheel, source distribution, and synthetic fixtures.

## Frozen recommended release

Use the official content release **2026-06-03**, not an unversioned “latest”
URL. The later `2026-06-11-ZENODO` maintenance release states that its database
content is unchanged from 2026-06-03 and attaches the same build.

| Field | Value |
|---|---|
| Official project | [antigenomics/vdjdb-db](https://github.com/antigenomics/vdjdb-db) |
| Release page | [2026-06-03](https://github.com/antigenomics/vdjdb-db/releases/tag/2026-06-03) |
| Exact asset | `https://github.com/antigenomics/vdjdb-db/releases/download/2026-06-03/vdjdb-2026-06-03.zip` |
| Release identifier | `2026-06-03` |
| Asset type | ZIP archive containing released tab-delimited VDJdb database tables, including the single-chain `vdjdb.txt` table |
| Published size | 41,017,409 bytes |
| Observed ZIP SHA256 | `9ae772461670d97040a6f2d9c719455bddce21fa9f16ef09e7a0bd35439cfc1` |
| Extracted `vdjdb.txt` size | 239,910,782 bytes |
| Extracted SHA256 | `3f5823f01d954751567a8d40870a284e384625242bce870b6390c3222a46568a` |
| Validated rows / unique CDR3 | 284,546 / 178,974 |
| Chain rows | TRA 122,243; TRB 162,303 |
| Citation | Goncharov M. et al. *VDJdb in the pandemic era: a compendium of T cell receptors specific for SARS-CoV-2*. Nature Methods (2022), DOI [10.1038/s41592-022-01578-0](https://doi.org/10.1038/s41592-022-01578-0) |
| Access/license | Open access; the upstream Zenodo metadata declares `AGPL-3.0-only` |

The database must remain an external runtime input. Its license and underlying
publication provenance should be cited, but the database must not be bundled in
the scRFU wheel, source distribution, repository, test fixtures, or manuscript
source-data archive.

## Executed evidence status

All 24 prespecified combinations across full Wells, GSE190905 and GSE157007
completed for CDR3/CDR3+V, nearest/threshold-qualified and
fractional/exclude-ambiguous policies. Four 1,000-permutation null models were
run per combination where the matched subset supported the statistic; sparse
undefined cases are recorded as skipped rather than optimized away.

The consolidated external source tables contain 24 sensitivity rows, 432
grouping-comparison rows, 96 null-summary rows and 96 evidence-score rows. Their
current SHA256 values are respectively:

- `78d9c9832559b6227a770d833cfc4e129a8fc80a2f4efd93ffed10b9fb34e5b2`
- `66ae878e3cf5a4825a27a77e449a2a26d4099529ea21701d708ed68704cedae0`
- `98aaf1994afca50330c3f917ea1a5d1245faaff7af1cd68add9f5d37a0842b65`
- `511b76dd8bf0968237c8ebb257dbe5e45419eef59163efc4b35b9abe68c95c27`

These results test coherence of external antigen annotations among distinct
sequences in the same RFU. They do not establish antigen specificity.

## Reproducible acquisition and configuration command

Replace only the external destination root. This command downloads the fixed
release, records its SHA256, extracts it, resolves exactly one `vdjdb.txt`, and
exports the variables expected by scRFU.

```bash
export VDJDB_RELEASE=2026-06-03
export SCRFU_VDJDB_ROOT=/path/to/external/vdjdb/2026-06-03
mkdir -p "$SCRFU_VDJDB_ROOT"
curl -fL \
  https://github.com/antigenomics/vdjdb-db/releases/download/2026-06-03/vdjdb-2026-06-03.zip \
  -o "$SCRFU_VDJDB_ROOT/vdjdb-2026-06-03.zip"
test "$(stat -c %s "$SCRFU_VDJDB_ROOT/vdjdb-2026-06-03.zip")" -eq 41017409
sha256sum "$SCRFU_VDJDB_ROOT/vdjdb-2026-06-03.zip" \
  > "$SCRFU_VDJDB_ROOT/vdjdb-2026-06-03.zip.sha256"
unzip -q "$SCRFU_VDJDB_ROOT/vdjdb-2026-06-03.zip" -d "$SCRFU_VDJDB_ROOT/release"
export VDJDB_PATH="$(find "$SCRFU_VDJDB_ROOT/release" -type f -name vdjdb.txt -print -quit)"
test -n "$VDJDB_PATH"
test "$(find "$SCRFU_VDJDB_ROOT/release" -type f -name vdjdb.txt | wc -l)" -eq 1
printf 'VDJDB_RELEASE=%s\nVDJDB_PATH=%s\n' "$VDJDB_RELEASE" "$VDJDB_PATH"
```

The command is retained for independent reproduction. Runtime locations must be
supplied by the user and are never embedded in tracked manifests.
