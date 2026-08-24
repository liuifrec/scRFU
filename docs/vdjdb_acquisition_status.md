# VDJdb acquisition status

Checked on 2026-08-24. `VDJDB_PATH` and `VDJDB_RELEASE` are unset, so no VDJdb
analysis was run and no database was downloaded automatically.

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
| Citation | Goncharov M. et al. *VDJdb in the pandemic era: a compendium of T cell receptors specific for SARS-CoV-2*. Nature Methods (2022), DOI [10.1038/s41592-022-01578-0](https://doi.org/10.1038/s41592-022-01578-0) |
| Access/license | Open access; the upstream Zenodo metadata declares `AGPL-3.0-only` |

The database must remain an external runtime input. Its license and underlying
publication provenance should be cited, but the database must not be bundled in
the scRFU wheel, source distribution, repository, test fixtures, or manuscript
source-data archive.

## Exact acquisition and configuration command

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

After those exports are present, rerun the final-scientific-validation session.
The reference validator will record the extracted table's SHA256, dimensions,
TRB counts, antigen fields, confidence/evidence fields, and duplicate evidence
before any matching or permutation analysis.
