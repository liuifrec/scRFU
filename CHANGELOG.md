# Changelog

## Unreleased

- Correct headerless upstream RFU input handling.
- Deduplicate exact CDR3 queries while reconstructing every input row stably.
- Retain ineligible and unassigned rows with explicit status fields.
- Add restartable, validated serial RFU chunking and cache reuse provenance.
- Add canonical receptor tables plus Wells, scirpy/AIRR, generic DataFrame, and
  Cell Ranger VDJ adapters.
- Add selective generic H5AD readers and schema-2 receptor caches.
- Reduce Wells receptor preparation memory use without loading expression.
- Add conventional repertoire metrics, refined RFU convergence metrics,
  sample-level RFU pseudobulk, RFU overlap, and phenotype coupling.
- Add dataset-independent plotting, sequence-matrix, downstream workflow, and
  benchmark-summary foundations.
- Add a local, version-labelled VDJdb reference contract with optional SHA256
  verification and explicit source-column mappings.
- Add exact CDR3, CDR3-plus-V, and genuinely paired evidence matching as a
  validated long table with stable evidence tiers.
- Add ambiguity-preserving sequence/row summaries, RFU antigen-coherence and
  global association metrics, size-preserving stratified nulls, simple receptor
  grouping baselines, biological-metadata recurrence summaries, and generic
  antigen-evidence plots.
- Add an offline VDJdb evidence workflow; no database content, live API, or
  silent download is included.

The recommended release version is 0.2.0. The package version is not changed
until repository release policy and placeholder project metadata are confirmed.
