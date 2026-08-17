# Methods-paper public API freeze

This document freezes the intended methods-paper interface; it does not declare
package version 1.0. Existing public functions are not removed in this phase.

## Stable for the methods paper

| Namespace | Symbols | Contract |
|---|---|---|
| `scrfu.pp` | `canonicalize_receptor_table`, `validate_receptor_table`, `receptor_schema`, normalization helpers | DataFrame in; canonical DataFrame or QC mapping out. Canonical receptor rows retain explicit row, cell, chain, sequence, gene, and source identifiers. |
| `scrfu.adapters` | `prepare_receptors`, `get_receptor_adapter`, `list_receptor_adapters`, named adapter functions, `AdapterResult` | User data in; canonical receptors, separate cell metadata, QC, adapter identity, and provenance out. MuData-like routing requires explicit modalities. |
| `scrfu.io` | selective H5AD readers, receptor-cache read/write/validate/migrate, `file_sha256`, `export_rfu_matrix` | Expression is not materialized by selective readers. Cache schema and checksums are validated. |
| `scrfu.tl` | `call_rfu_table`, `call_rfu`, `call_rfu_repo`, `validate_airr`, `rfu_metrics`, `repertoire_metrics`, `rfu_pseudobulk`, `rfu_overlap`, `rfu_phenotype_coupling`, `rfu_sequence_matrix` | Table-level RFU execution is canonical; AnnData functions delegate. Assignment and weighting policies are explicit. |
| `scrfu.pl` | existing RFU summary plots and generic downstream/antigen plots | `ax=None` accepted, Matplotlib `Axes` returned, no `show()`. Plotting is optional. |

The offline VDJdb reference, exact matching, evidence summaries, coherence,
grouping comparison, and permutation functions exposed through `scrfu.tl` are
also stable as evidence-analysis contracts. A match is annotation evidence, not
proof of specificity.

## Experimental during Month 1

The following are public and tested synthetically, but remain experimental until
real longitudinal or multi-cohort validation freezes their scientific defaults:

- `reference_coverage`;
- `validate_longitudinal_design`, `rfu_longitudinal_matrix`,
  `longitudinal_similarity`, `summarize_longitudinal_similarity`,
  `donor_retrieval`, `rfu_longitudinal_dynamics`,
  `longitudinal_compartment_comparison`, `bootstrap_longitudinal_statistic`, and
  `permute_longitudinal_labels`;
- `FrozenRFUReference`, `validate_frozen_reference`, `transfer_cohort`,
  `harmonize_cohort_metadata`, and `create_heldout_validation_manifest`;
- deterministic subsampling, abundance resampling, stability benchmarks,
  threshold sensitivity, and comparator registration/representations;
- `max_workers` and `executor` on RFU calls (serial remains the default).

Experimental result dataclasses and their current field names are documented and
tested, but may gain fields before the release candidate. Fields will not be
silently reinterpreted.

Public entry points currently use `TypeError`, `ValueError`, `KeyError`,
`FileNotFoundError`, and actionable `ImportError` for input/configuration
problems. Backend-specific capability and execution exceptions remain available
from their established backend modules but are not promoted as new top-level
symbols in Month 1. CLI commands must express the same defaults and validation
as their Python counterparts; new longitudinal/transfer workflows are examples,
not yet core CLI subcommands.

## Compatibility aliases

`rfu_summary`, `aggregate_rfu`, adapter aliases (`wells`, `cellranger`,
`tenx_vdj`, `airr`, `dataframe`), legacy RFU output columns (`rfu_id`,
`rfu_label`, `pass_thr`), and AnnData attachment names remain supported.
`rfu_id_nearest`, `rfu_label_nearest`, and `rfu_pass_threshold` make their
semantics explicit. No removal warning is added in Month 1.

## Internal only

Underscore-prefixed helpers, `analysis_frame`, low-level chunk runners,
reconstruction helpers, VDJdb matching internals, visual alignment helpers, and
module implementation paths outside the public namespaces are internal. Users
must not depend on chunk-directory internals except for the documented manifest
schema.

## Future or deferred

Paired-chain functional assignment, a BCR functional reference, inferential
longitudinal models, automatic cohort ontology inference, and expression-level
batch integration are deferred. MuData input routing does not imply those
capabilities.

## Frozen terminology

- **receptor row**: one selected receptor record; not necessarily one cell;
- **cell**: a unique biological cell identifier;
- **sequence**: an exact amino-acid CDR3 unless another sequence field is named;
- **clonotype**: only a supplied clonotype identity, never silently synthesized
  from CDR3;
- **functional unit**: a label from the frozen RFU reference;
- **sample**: the independent biological specimen represented by one matrix row;
- **phenotype**: a descriptive cell/sample label;
- **compartment** and **timepoint**: explicit repeated-measures design fields;
- **nearest assignment**: the best returned reference unit irrespective of
  threshold;
- **threshold-qualified assignment**: a nearest assignment meeting the frozen
  similarity threshold;
- **reference coverage**: descriptive eligibility/assignment/threshold counts,
  not a probabilistic confidence interval or validated OOD classifier.

## Symbol-by-symbol export classification

All names below are current `__all__` exports. “Stable” refers to the
methods-paper compatibility intent; “experimental” requires real validation as
described above.

| Namespace | Stable exports | Experimental exports |
|---|---|---|
| `scrfu.pp` | `OPTIONAL_RECEPTOR_COLUMNS`, `RECEPTOR_SCHEMA_VERSION`, `RECOGNIZED_CHAINS`, `REQUIRED_RECEPTOR_COLUMNS`, `canonicalize_receptor_table`, `normalize_chain`, `normalize_productive`, `receptor_schema`, `validate_receptor_table` | — |
| `scrfu.adapters` | `ADAPTER_API_VERSION`, `AdapterResult`, `ReceptorAdapter`, `adapt_airr_dataframe`, `adapt_cellranger_vdj`, `adapt_wells_tcr_ir`, `get_receptor_adapter`, `list_receptor_adapters`, `prepare_receptors` | MuData-specific `modality`/`metadata_modality` arguments only |
| `scrfu.io` | `H5ADDataFrameInfo`, `RECEPTOR_CACHE_SCHEMA_VERSION`, `ReceptorCacheData`, `UnsupportedH5ADLayout`, `export_rfu_matrix`, `file_sha256`, `migrate_wells_receptor_cache`, `read_h5ad`, `read_h5ad_dataframe`, `read_h5ad_obs`, `read_h5ad_shape`, `read_receptor_cache`, `source_fingerprint`, `validate_receptor_cache`, `write_h5ad`, `write_receptor_cache` | — |
| `scrfu.tl` | `AntigenContextResult`, `AntigenPermutationResult`, `RFUOverlapResult`, `RFUPseudobulkResult`, `RFUTableResult`, `VDJdbEvidenceSummary`, `VDJdbReference`, `aggregate_rfu`, `annotate_vdjdb`, `call_rfu`, `call_rfu_repo`, `call_rfu_table`, `compare_antigen_groupings`, `global_antigen_coherence`, `load_vdjdb_reference`, `normalize_vdjdb_cdr3`, `normalize_vdjdb_v_gene`, `repertoire_metrics`, `rfu_antigen_abundance`, `rfu_antigen_coherence`, `rfu_antigen_permutation_test`, `rfu_metrics`, `rfu_overlap`, `rfu_phenotype_coupling`, `rfu_pseudobulk`, `rfu_sequence_matrix`, `rfu_summary`, `summarize_antigen_context`, `summarize_vdjdb_evidence`, `validate_airr`, `validate_vdjdb_reference` | `CohortHarmonizationResult`, `ComparatorRepresentation`, `FrozenRFUReference`, `HeldOutValidationManifest`, `LongitudinalCompartmentResult`, `LongitudinalDesign`, `LongitudinalDynamicsResult`, `LongitudinalResamplingResult`, `RFULongitudinalResult`, `StabilityBenchmarkResult`, `TransferCohortResult`, `benchmark_representation_stability`, `bootstrap_longitudinal_statistic`, `create_heldout_validation_manifest`, `deterministic_subsample`, `donor_leave_one_out`, `donor_retrieval`, `harmonize_cohort_metadata`, `list_comparators`, `longitudinal_compartment_comparison`, `longitudinal_similarity`, `multinomial_abundance_resample`, `permute_longitudinal_labels`, `reference_coverage`, `register_comparator`, `repertoire_representation`, `rfu_longitudinal_dynamics`, `rfu_longitudinal_matrix`, `shuffle_input_order`, `summarize_longitudinal_similarity`, `threshold_sensitivity`, `transfer_cohort`, `validate_frozen_reference`, `validate_longitudinal_design` |
| `scrfu.pl` | `antigen_permutation_distribution`, `repertoire_metric_comparison`, `rfu_antigen_bubble`, `rfu_antigen_coherence`, `rfu_antigen_heatmap`, `rfu_bar`, `rfu_convergence`, `rfu_heatmap`, `rfu_metric_heatmap`, `rfu_overlap_heatmap`, `rfu_phenotype_heatmap`, `rfu_score_hist` | — |

There are no Month 1 public exports classified for removal. Underscore-prefixed
symbols and unexported implementation classes remain internal.
