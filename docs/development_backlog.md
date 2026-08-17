# Issue-ready development backlog

Each item lists priority; rationale; dependencies; implementation; test;
acceptance; manuscript support; and status.

## 1. Core assignment

- **P0 — Frozen original-RFU parity table.** Rationale: technical fidelity is
  non-negotiable. Depends on version-pinned external artifacts. Implement the
  comparison workflow and mismatch export; test synthetic mismatches; accept
  zero unexplained ID/threshold differences and tolerance-bounded scores.
  Supports Figure 1 parity. **Blocked by artifact run.**
- **P1 — Real parallel scaling.** Depends on parity and bounded caches. Record
  workers/chunks/runtime/RSS; test manifest reconciliation; accept identical
  scientific output and no corrupted completed chunks. Figure 1 performance.
  **Implementation unblocked; run blocked.**

## 2. Data integration

- **P0 — Installed-MuData compatibility fixture.** Rationale: validate the
  optional interface, not just duck typing. Add optional dependency CI test;
  accept untouched modalities, aligned cells, actionable failures. Figure 1
  inputs. **Unblocked.**
- **P0 — Model/chain compatibility guard.** Depends on frozen reference schema.
  Reject TCR-reference/BCR misuse; test every chain/model combination; accept no
  silent cross-receptor assignment. Figure 1 QC. **Unblocked.**

## 3. Longitudinal analysis

- **P0 — Prespecify trajectory thresholds and sensitivity.** Depends on cohort
  design review. Export all trajectories/missingness; test boundaries; accept
  stable labels under declared settings. Figure 2. **Blocked by scientific
  review.**
- **P0 — Leave-one-timepoint donor retrieval benchmark.** Depends on public and
  deep longitudinal inputs; reuse identical comparator interface; accept frozen
  metrics with candidate counts and donor-level uncertainty. Figure 2.
  **Blocked by datasets.**

## 4. Cross-cohort transfer

- **P0 — Register frozen reference and held-out cohort.** Depends on artifact
  hashes and dataset verification. Save manifest before evaluation; test role
  overlap; accept immutable ID and exact data hashes. Figure 3. **Blocked.**
- **P0 — Cross-cohort coverage/harmonization report.** Explicit maps only;
  accept all unmapped values and missingness reported. Figure 3. **Blocked by
  cohort selection.**

## 5. BCR architecture

- **P2 — BCR go/no-go review.** Requires dedicated reference, public cohort,
  isotype/SHM/family semantics, fidelity, and held-out plan. No implementation
  until approved. Optional Figure 4. **Blocked by design gate.**

## 6. Benchmarks

- **P0 — 1k/10k/25k scaling and cache reuse.** Depends on public cache and
  resource window. Test manifest schema; accept complete runtime/RSS/hash
  tables. Figure 1. **Run blocked.**
- **P1 — Robustness perturbation suite expansion.** Add clone-frequency and
  source-count-aware depth modes only when semantics are defined; deterministic
  tests; accept source units correctly labelled. Figure 1/supplement.
  **Partially unblocked.**

## 7. Visualization

- **P1 — Longitudinal and transfer plotting layer.** Depends on frozen source
  table schemas. Matplotlib only; test axes/empty/top-N; accept no implicit
  inferential annotations. Figures 2/3. **Blocked by schema review.**

## 8. Documentation and release

- **P0 — Resolve package/tag version history.** Maintainer decision required;
  metadata/artifact tests; accept package, tag, changelog, and release agreement.
  Submission gate. **Blocked.**
- **P1 — Clean-wheel longitudinal tutorial.** Depends on API freeze; offline
  synthetic workflow; accept Python 3.10–3.12. Supplement. **Unblocked.**

## 9. Manuscript figures

- **P0 — Panel-to-source-table index.** Depends on final analyses; test file/hash
  completeness; accept one immutable source table per panel. Figures 1–3.
  **Blocked by results.**

## 10. Public validation datasets

- **P0 — Verify two TCR cohorts and one held-out role.** Rationale:
  transferability cannot rest on one atlas. Record fields, license, confounding,
  preprocessing, and role before analysis; accept independent public evidence.
  Figures 2/3. **Blocked by dataset verification.**
- **P1 — Decide small public fixture.** Depends on license/redistribution review;
  accept tiny documented download or synthetic fallback. Release/tutorial.
  **Blocked by licensing review.**
