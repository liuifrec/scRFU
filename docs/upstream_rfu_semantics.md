# Upstream RFU Semantics and Backend Capabilities

## Scope and inspected implementations

This note records the behavior that scRFU must preserve before adding
large-atlas orchestration. It distinguishes two implementations:

- The official public [s175573/RFU repository](https://github.com/s175573/RFU)
  at commit `ff7a1ea6aca28444a7ff8ec7ef09745e7341f99d`. Its public `RFU.R` defines
  `EncodeRepertoire()` and `AssignRFUs()`, but not the map-aware functions.
- A local development extension of `RFU.R` inspected for behavioral comparison.
  It additionally defines `AssignRFUs_with_map()` and
  `RFUbatch_with_maps()`. Those functions are not treated as public upstream
  capabilities and are not required by standard scRFU operation.

The local development checkout is configuration, not package data. No local
absolute path belongs in source, tests, package metadata, examples, or public
documentation. Examples must use a neutral path such as `/path/to/RFU`.

## Official public behavior

### `getTrimerPos()` and `EncodeRepertoire()`

`EncodeRepertoire(ff)` reads a tab-separated file with `header=FALSE` and:

1. retains input rows whose first column starts with `C`;
2. obtains position-labelled amino-acid trimers with `getTrimerPos()`, using
   `st=2` and `ed=4`;
3. removes the position labels and retains trimers present in the rows of the
   loaded `fit` object;
4. uses the corresponding `fit` row for one retained trimer or the column-wise
   mean embedding for multiple retained trimers; and
5. returns one embedding row per retained input row.

The active implementation reads the second input column but does not use its
TRBV value in encoding or assignment. A historical, commented implementation
did include V-gene CDR embeddings; it is not the behavior of public
`AssignRFUs()`.

Input order and duplicate sequences are preserved after the first-column
`^C` filter. Consequently, repeated CDR3 sequences contribute repeated rows and
repeated counts. The function does not expose stable input row identifiers.

### `AssignRFUs()`

For each encoded repertoire row, `AssignRFUs()` computes Spearman correlations
against all centroid rows in `CL$centers`. It returns:

- `TCR`: the one-based index of the first maximum-correlation centroid for each
  encoded row;
- `COR`: that row's maximum correlation;
- `N`: the number of encoded rows with maximum correlation below `THR`; and
- `RFU`: counts of all `TCR` assignments, divided by the number of encoded rows
  and multiplied by 10,000. With `normalize=TRUE`, each value is additionally
  divided by the corresponding reference cluster size from `table(CL$cl)`.

The threshold does **not** remove low-correlation rows from `TCR` or from the
`RFU` abundance vector. It only determines `N` (and can be used by a wrapper to
derive a pass flag). Any scRFU implementation that changes this rule would not
be scientifically equivalent to the official function.

### Reproducible public capability

An official public checkout containing `RFU.R`,
`trimerMDSfit_small.Rdata`, and `km5000noMax.Rdata` is sufficient to reproduce:

- the standard `AssignRFUs()` RFU abundance vector;
- the per-encoded-row winning RFU index and maximum correlation; and
- the official threshold-derived miss count.

Stable per-cell reconstruction is implemented in scRFU without private R code.
The Python layer assigns `input_row_id`, applies the same `^C` eligibility rule,
queries by `unique_sequence_id`, and reconstructs with a validated many-to-one
mapping. It never joins unrestricted assignment results on CDR3 text.

The previous wrapper wrote a column header named `CDR3` into the intermediate
file. Because `EncodeRepertoire()` reads with `header=FALSE` and `CDR3` itself
matches `^C`, that header became a synthetic encoded repertoire row. The current
wrapper writes the upstream intermediate without column names.

## Optional local map-aware extension

`AssignRFUs_with_map()` repeats the official correlation, assignment, counting,
normalization, and threshold calculations. It also reads the first three input
columns as CDR3, TRBV, and optional frequency and returns `MAP` with `cdr3_aa`,
`trbv`, `freq`, `rfu`, `max_cor`, and `pass_thr`.

The extension deliberately does not apply the `^C` filter to its metadata
table. It pairs the first `min(number of metadata rows, number of encoded rows)`
rows by position. Therefore its map becomes shifted when `EncodeRepertoire()`
drops an earlier non-`C` row. Truncation prevents a length error but does not
guarantee stable identity. The map is scientifically equivalent only when the
metadata rows and encoded rows are already in one-to-one order.

`RFUbatch_with_maps()` enumerates matching TSV/TXT files, calls
`AssignRFUs_with_map()` for each, column-binds the standard RFU vectors, writes
one mapping file per sample, and writes `RFU_matrix.tsv`. It provides batch I/O,
not a different assignment algorithm, and inherits the mapping limitation.

These functions may be used only as optional capabilities when a user supplies
a checkout that contains them and their license permits use. Public scRFU must
not require them.

In `map_aware` mode, scRFU does not accept the extension's `MAP` row identifiers
as authoritative. It submits only eligible queries with explicit
`unique_sequence_id` values and reconstructs from the length-checked `TCR` and
`COR` vectors. A focused integration test places an ineligible row before two
eligible rows and verifies that both later assignments retain their original
identifiers.

## Backend resolution and capability design

The implemented backend behavior is:

1. Resolve the checkout from an explicit `rfu_dir` argument when supplied.
2. Otherwise resolve it from the `RFU_DIR` environment variable.
3. Otherwise raise an actionable error that describes both supported options
   and uses `/path/to/RFU` as its example.
4. Validate the three files required for standard mode.
5. Read `RFU.R` without executing it and detect exact top-level function
   definitions for:
   `AssignRFUs`, `AssignRFUs_with_map`, and `RFUbatch_with_maps`.
6. Select `standard` mode by default. Standard mode requires only
   `AssignRFUs`. Optional `map_aware` mode requires
   `AssignRFUs_with_map`; requesting unavailable map-aware or batch-map behavior
   raises a capability error that lists the detected functions and explains
   that the official public checkout supports standard mode.

Capability tests construct temporary fake RFU directories and synthetic
`RFU.R` files. They must not depend on a developer checkout or an absolute
machine path.

Each run's provenance should record:

- the resolved backend path;
- SHA-256 hashes of `RFU.R`, `trimerMDSfit_small.Rdata`, and
  `km5000noMax.Rdata`;
- detected RFU functions;
- selected backend mode;
- the scRFU version; and
- the wrapper path/hash, R executable, timestamp, and run parameters already
  tracked by scRFU.

`standard` is the canonical default even when optional functions are detected;
capability presence never silently changes the selected scientific mode.

Public-standard integration tests use only `Rscript`, `RFU_DIR`, the three
required backend files, and `AssignRFUs()`. Optional map-aware tests use the
separate `SCRFU_MAP_AWARE_RFU_DIR`; they cannot make a standard test available
or unavailable. Run the two groups independently with neutral user paths:

```bash
RFU_DIR=/path/to/RFU python -m pytest -q \
  tests/test_integration_rfu_repo.py -k "not map_aware"

SCRFU_MAP_AWARE_RFU_DIR=/path/to/map-aware-RFU python -m pytest -q \
  tests/test_integration_rfu_repo.py -k map_aware
```

Skip messages distinguish an unavailable `Rscript`, an unset or incomplete
environment-specific checkout, and a missing requested capability.

## Stable standard mapping and equivalence

The standard mapping path preserves every extracted row and its original order.
Rows failing `^C` remain in the reconstructed result with
`eligibility_status="ineligible_cdr3_not_starting_c"` and no RFU assignment.
Eligible rows receive a `unique_sequence_id`. With deduplication enabled, the
first occurrence of each exact `cdr3aa` value is submitted and all occurrences
map back through that identifier. TRBV is preserved per original row but is not
part of the candidate deduplication key because active `EncodeRepertoire()` does
not use it.

The external-backend integration test compares repeated-row calls, unique-CDR3
calls, and reconstructed repeated rows. For every eligible row it requires the
same nearest RFU ID and label, identical threshold status, and RFU score within
absolute tolerance `1e-12`. The test includes one CDR3 repeated with different
TRBV calls, unique sequences, a non-`C` sequence, and at least two below-threshold
sequences. It separately verifies sequence-level equivalence, restored
repertoire multiplicity, unique-query upstream `N`, and multiplicity-weighted
reconstructed miss count.

## Restartable unique-sequence chunking

The execution order is fixed: original rows receive `input_row_id`; the `^C`
eligibility rule is applied; eligible rows are deduplicated by exact `cdr3aa`;
the ordered unique-sequence table is divided into deterministic chunks; chunk
outputs are validated; unique-sequence results are concatenated in chunk order;
and the existing validated many-to-one map restores original order and
multiplicity. Original per-cell rows are never chunked before equivalent CDR3
queries are deduplicated.

`chunk_size=None` keeps the existing one-call path. A positive `chunk_size`
enables serial restartable execution. `resume=True` is the default;
`resume=False` reruns all chunks without reusing cached output, and
`force_recompute=True` takes precedence and forces recomputation. The Wells
workflow chooses 20,000 unique CDR3s by default without changing small calls to
the core API.

### Deterministic identity

The run ID is a SHA-256 digest over canonical JSON containing ordered
`unique_sequence_id` values, ordered CDR3 sequences, backend mode, threshold,
deduplication key, eligibility rule, hashes of `RFU.R`, both atlas files, and
the R wrapper, the wrapper schema and scRFU versions, chunk size, and extra R
arguments. Runtime timestamps, work directories, and temporary names are not
included.

For each chunk, scRFU records its zero-based index, start and exclusive end
offsets, expected row count, ordered identifiers and their hash, and the hash of
the exact TSV input bytes. Its ID is
`<first-12-run-hash>-<five-digit-index>-<first-12-input-hash>`.

### Run directory and cache validation

Chunked runs use this stable structure:

```text
<workdir>/runs/<run_id>/
  run_manifest.json
  chunks/
    chunk_00000/
      input.tsv
      output.tsv
      stdout.log
      stderr.log
      manifest.json
```

Manifest and output writes use temporary sibling files followed by atomic
replacement. A chunk is a cache hit only when the manifest parses, uses the
supported schema, is complete, and matches the run/chunk IDs, index, offsets,
expected count, input and ordered-identifier hashes, wrapper and backend hashes,
mode, threshold, eligibility rule, and deduplication key. The output must exist,
match its recorded hash, contain every required column, contain the expected
number of rows, have unique identifiers, and match the exact expected identifier
order. The presence of `output.tsv` alone is never sufficient.

Invalid entries are retained as diagnostics, classified with an invalidation
reason, and recomputed. A failed chunk retains stdout, stderr, its exit code,
and a failed manifest; earlier valid chunks remain intact. A resumed run can
therefore reuse earlier chunks and restart at the failed chunk. Concatenation
occurs only after every chunk validates, in chunk-index order, with a final
global identifier/order check.

Synthetic tests compare unchunked execution with chunk sizes 1, 2, and a
non-divisible final chunk and require identical IDs, labels, statuses,
identifiers, and row order, with RFU scores equal to absolute tolerance
`1e-12`. The official-compatible integration test independently compares
chunked and unchunked standard `AssignRFUs()` execution under the same
tolerance. These checks demonstrate sequence-assignment and reconstructed-row
equivalence for the tested inputs; they do not assert behavior for a different
upstream implementation.

## Wells atlas `TCR_IR` adapter

The Wells atlas processing workflow used Scirpy 0.10.1, prefixed TCR fields with
`TCR-`, and collected them into a flattened `TCR_IR` table. Its public processing
code stores that table in `adata.obsm["TCR_IR"]`; some distributed objects store
the same table in `adata.uns["TCR_IR"]`.

scRFU accepts both locations when `airr_key="TCR_IR"`, with `uns` taking
precedence when both exist. Recognized primary VDJ fields are:

- CDR3 amino acid: `TCR-IR_VDJ_1_junction_aa` or
  `IR_VDJ_1_junction_aa`, with legacy `..._cdr3` aliases;
- V gene: `TCR-IR_VDJ_1_v_call` or `IR_VDJ_1_v_call`;
- optional locus: `TCR-IR_VDJ_1_locus` or `IR_VDJ_1_locus`; and
- optional productive flag: `TCR-IR_VDJ_1_productive` or
  `IR_VDJ_1_productive`.

The adapter uses an explicit cell column when available. Otherwise the table
must have the same number and exact index order as `adata.obs_names`, or have a
default positional index with one row per observation. Ambiguous alignment and
duplicate explicit cell identifiers are errors. Missing values (including the
literal `"nan"` emitted by the atlas preprocessing) and non-target chains are
removed. Non-`C` CDR3 sequences remain in the feature table; the backend retains
and explicitly marks these rows instead of submitting them to
`EncodeRepertoire()`. The adapter returns the existing
`cell_id`/`cdr3aa`/`trbv` feature schema and does not mutate the AnnData object.
