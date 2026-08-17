# Methodological novelty boundary

## Original RFU contribution

The original RFU work defines the trimer representation, the RFU reference,
centroids or associated reference artifacts, and the original assignment
method. Exact parity with that implementation is a technical requirement for
scRFU, not a claim that scRFU invented RFUs. Users must cite the original RFU
method when using its assignments.

## scRFU contribution

scRFU contributes the transferable analysis and execution layer around a frozen
RFU reference:

- standardized AIRR, scirpy, AnnData, selective-H5AD, Cell Ranger, Wells, and
  bulk receptor-table inputs;
- canonical receptor rows with metadata kept in a separate, validated table;
- exact-sequence deduplication followed by stable cell/receptor-row
  reconstruction;
- restartable chunk manifests and deterministic optional chunk parallelism;
- hashes, thresholds, eligibility rules, assignment modes, and source-table
  provenance;
- explicit nearest versus threshold-qualified assignments and descriptive
  low-reference-similarity coverage diagnostics;
- sample-level RFU representations, repertoire baselines, phenotype linkage,
  and offline version-pinned VDJdb evidence analysis;
- repeated-measures design validation, longitudinal similarity, donor
  retrieval, explicit RFU trajectory classification, and donor-block
  resampling;
- frozen-reference transfer summaries and explicit metadata harmonization;
- deterministic robustness utilities and shared sample-representation
  comparator interfaces.

These are methods claims only after synthetic correctness, original-RFU parity,
bounded public validation, independent-cohort validation, and held-out testing
have each been demonstrated at the level stated in the manuscript.

## Future or gated contribution

Receptor-specific BCR functional units, inferential longitudinal models, and a
cross-receptor synthesis are not current capabilities. BCR remains gated on a
separate reference, technical fidelity, public-cohort analysis, interpretable
maturation/isotype behavior, and held-out validation. The canonical schema's
ability to store immunoglobulin chains is data-model support only.
