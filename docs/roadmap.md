# Roadmap

scRFU is being developed as a transferable frozen-reference receptor-state
framework, not solely as an RFU wrapper. The three-month execution criteria are
in [`three_month_plan.md`](three_month_plan.md).

## Month 1: method and API freeze — current

- freeze receptor, assignment, coverage, longitudinal, transfer, comparator,
  and provenance terminology;
- preserve original-RFU semantics and backwards compatibility;
- provide optional deterministic parallel chunks with serial default;
- validate repeated-measures and frozen-reference contracts synthetically;
- freeze benchmark datasets, metrics, figures, source tables, and BCR gate.

Real longitudinal cohorts and full public cohorts are deliberately not run in
this phase.

## Month 2: benchmarks and applications

- original-RFU parity and bounded 1k/10k/25k scaling;
- downsampling, order, threshold, chunk, and worker robustness;
- prespecified deep longitudinal methods demonstration;
- independent public TCR and paired single-cell applications;
- completely held-out transfer evaluation;
- BCR implementation only if every architecture gate passes.

## Month 3: manuscript and release

- reproducible figures and source tables;
- clean-install tutorials, public test material, CI, and API reference;
- approved versioned release and DOI archive;
- manuscript, methods, cover letter, and reproducibility report.

## Explicitly deferred

Inferential longitudinal models, automatic ontology mapping, expression batch
integration, cross-receptor synthesis, and BCR functional units without a
validated receptor-specific reference.
