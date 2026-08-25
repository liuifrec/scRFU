# Changelog

## Unreleased

Development version: `0.4.0rc1`. The intended next public release is 0.4.0;
no release or tag has been created.

### Added

- Canonical receptor rows; Wells, AIRR/scirpy, generic DataFrame, and Cell
  Ranger VDJ adapters; selective H5AD extraction; and portable checksummed
  receptor caches.
- Dataset-independent table RFU engine, exact-CDR3 deduplication, stable row
  reconstruction, restartable chunks, and optional deterministic process/thread
  chunk workers with serial default.
- Conventional repertoire and RFU convergence metrics, sample pseudobulk,
  overlap, phenotype coupling, generic plotting, sequence-frequency matrices,
  downstream workflow, and benchmark manifest summary.
- Local version-labelled VDJdb references with SHA256 validation, exact
  long-format evidence matching, ambiguity-preserving summaries, RFU antigen
  coherence, grouping baselines, stratified permutation nulls, and offline
  workflow/plots.
- Explicit reference-coverage diagnostics; repeated-measures design validation;
  longitudinal matrices, similarity, donor retrieval, trajectory classes,
  compartment comparisons, and donor-block bootstrap/permutation foundations.
- Immutable frozen-reference and held-out-validation manifests, target-cohort
  transfer summaries, explicit metadata harmonization, deterministic robustness
  utilities, and a common comparator representation interface.
- Optional MuData-like modality routing without a mandatory MuData dependency.
- Experimental BCR canonicalization, heavy/light pairing, QC, conservative
  receptor-state features, and synthetic AnnData/MuData-like integration tests;
  no BCR functional-unit reference is defined.
- A redistributable synthetic end-to-end tutorial fixture, generated scaling
  benchmark, and advisory performance-regression comparison utility.

### Changed

- Official `AssignRFUs()` standard mode remains canonical; AnnData APIs delegate
  to the table engine.
- Metadata remain separate from receptor rows, and nearest assignments are
  distinguished from threshold-qualified assignments.
- Chunk provenance now records worker count and executor. Completed chunks are
  reusable across scientifically equivalent serial and parallel orchestration.
- Plotting imports remain lazy/optional for non-plotting core use.
- Project metadata now identifies the public repository and maintainer instead
  of placeholders.

### Fixed

- Headerless upstream RFU input no longer treats a CDR3 header as a synthetic
  sequence.
- Ineligible and upstream-unassigned rows are retained explicitly.
- Unsafe many-to-many CDR3 merges and unstable sequence-to-cell reconstruction
  were replaced by identifier-validated mappings.
- VDJdb CDR3+V evidence queries now remain distinct from CDR3-based RFU
  sequence identity, including heterogeneous or missing V calls without an
  arbitrary V-gene collapse.
- Missing R executables and malformed RFU wrapper outputs now produce explicit,
  actionable failures rather than ambiguous subprocess/parser errors.
- Map-aware identifier shifting, location-dependent subprocess tests, Python
  3.10 datetime handling, and large-H5AD expression materialization were
  corrected.

### Deprecated

- No public symbol is removed in this phase. Legacy `rfu_id`, `rfu_label`, and
  `pass_thr` remain compatibility fields alongside explicit nearest/threshold
  aliases. Any future removal will use `FutureWarning`, name the replacement,
  and document a target version.

### Migration notes

- Regenerate outputs produced by the older header-writing RFU wrapper.
- Treat nearest assignment and threshold-pass status as distinct; threshold
  failure is a low-reference-similarity flag, not calibrated confidence or OOD.
- Expect ineligible and unassigned receptor rows to remain in per-row outputs.
- Legacy Wells caches can be migrated to the generic receptor-cache schema.
- Serial execution remains the default. Parallelism requires explicit chunking
  and `max_workers > 1`.
- The development version is `0.4.0rc1` because historical `v0.2.0` and
  `v0.3.0` tags already exist. Do not infer a release solely from this section.
