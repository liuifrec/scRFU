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

Stable per-cell reconstruction is possible in scRFU without private R code when
the wrapper writes a validated, uniquely identified input table, applies the
same `^C` eligibility rule before calling R, and reconstructs output by retained
row position rather than by CDR3 text. Scientific-equivalence tests must cover
duplicates, filtered sequences, ties, thresholds, and empty or malformed input
before that reconstruction replaces existing behavior.

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

## Backend resolution and capability design

The backend design for the next implementation stage is:

1. Resolve the checkout from an explicit `rfu_dir` argument when supplied.
2. Otherwise resolve it from the `RFU_DIR` environment variable.
3. Otherwise raise an actionable error that describes both supported options
   and uses `/path/to/RFU` as its example.
4. Validate the three files required for standard mode.
5. Source `RFU.R` in an isolated R process and detect functions with exact
   function checks. Record at least:
   `AssignRFUs`, `AssignRFUs_with_map`, and `RFUbatch_with_maps`.
6. Select `standard` mode by default. Standard mode requires only
   `AssignRFUs`. Optional `map_aware` mode requires
   `AssignRFUs_with_map`; requesting unavailable map-aware or batch-map behavior
   raises a capability error that lists the detected functions and explains
   that the official public checkout supports standard mode.

Capability tests must construct temporary fake RFU directories and synthetic
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

Capability detection, chunk orchestration, and result reconstruction are design
requirements documented here, not changes included in the Wells-adapter stage.
The public RFU result schema remains unchanged during this stage.

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
literal `"nan"` emitted by the atlas preprocessing), non-target chains, and CDR3
sequences that official `EncodeRepertoire()` would discard are removed before
backend execution. The adapter returns the existing
`cell_id`/`cdr3aa`/`trbv` feature schema and does not mutate the AnnData object.
