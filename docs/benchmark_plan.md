# Benchmark plan

## Technical fidelity

Compare scRFU against frozen original-RFU artifacts at row level for RFU ID,
score, nearest assignment, threshold status, duplicates, non-C/ineligible
sequences, multiple-chain rejection/selection, input order, chunk size, and
official versus separately configured optional backend mode. Record every
artifact SHA256 and an explicit mismatch tolerance (zero for identifiers and
threshold booleans; documented numerical tolerance for scores).

## Scaling and execution

Planned public/cache-based sizes are 1,000, 10,000, and 25,000 source cells plus
an explicit user size. Record source cells, receptor rows, eligible rows, unique
sequences, chunks, workers, deduplication ratio, wall and CPU time, peak RSS,
cache/output sizes, cache-reuse time, and serial/parallel efficiency. The full
atlas is not an automatic benchmark.

## Robustness perturbations

- whole-cell, exact-sequence, sample, and donor leave-one-out subsampling;
- abundance-vector multinomial resampling (never labelled physical read
  downsampling without source read counts);
- threshold sensitivity, input-order shuffle, and chunk-size variation;
- later: clone-frequency and source-count-aware sequencing-depth downsampling.

Report RFU abundance/rank/richness, convergence, phenotype coupling, reference
coverage, donor retrieval, and longitudinal similarity using Spearman, cosine,
Jaccard, top-k overlap, absolute error, and rank displacement where defined.
Every random operation stores its seed and selected-unit provenance.

## Acceptance

Serial, parallel, chunk-size, and processing-order outputs must be scientifically
identical. Undefined metrics stay undefined. Scaling claims require repeated
runs or clearly labelled single-run measurements and complete environment
manifests.
