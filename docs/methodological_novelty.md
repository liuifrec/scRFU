# Methodological novelty boundary

## Original RFU contribution

The original RFU project defines the RFU representation/model, trimer
transformation, canonical assignment method, centroids/reference artifacts,
threshold semantics, and their biological motivation. scRFU reproduces this
assignment exactly; it does not claim to have invented RFUs or redistribute the
upstream assets. Users must cite and obtain the original RFU implementation
separately.

## Scirpy contribution

Scirpy supplies the AIRR/scverse receptor data model, chain indexing and QC,
general clonotype definitions, distance calculations, and broad immune-
repertoire ecosystem. scRFU uses those public conventions and compares against
genuine Scirpy clonotypes. It does not relabel Scirpy functionality as its own.

## scRFU contribution

scRFU contributes a reproducible method/execution layer for applying a frozen
RFU reference to single-cell and cohort receptor data:

- canonical receptor adapters and exact-CDR3 query deduplication with stable
  receptor-row reconstruction;
- deterministic restartable chunking, process parallelism, cache validation,
  artifact hashing, and explicit reference-coverage diagnostics;
- frozen-reference transfer without cohort-specific refitting;
- chain-aligned AnnData/MuData AIRR annotations that preserve every source
  chain, plus ambiguity-aware opt-in cell summaries;
- portable storage provenance, serialization/subsetting behavior, and guarded
  concatenation across compatible references;
- RFU-specific pseudobulk, convergence, overlap, phenotype coupling,
  longitudinal representation, robustness, and comparator analyses;
- version-pinned external VDJdb annotation/coherence analysis that keeps RFU
  CDR3 identity separate from strict chain+CDR3+V evidence-query identity.

## Method-level assessment

The remaining novelty is more than packaging: exact reconstruction,
large-repertoire execution/recovery, frozen cross-cohort feature mapping, and a
chain-aligned scverse receptor-state schema together define a transferable
analysis method that neither original RFU nor Scirpy alone provides. The
contribution is nonetheless conditional on the original RFU reference and
complementary to Scirpy. It should be presented as a scalable interoperable RFU
method, not as a new receptor similarity model or a replacement for clonotype
analysis.

BCR functional units, antigen specificity, deep longitudinal inference, and
supervised outcome prediction are outside this contribution.
