# Release-candidate scorecard

Candidate `0.4.0rc1`, assessed 2026-08-25. PASS means demonstrated for the
current local tree; it does not authorize a release. PARTIAL identifies an
honest boundary, not a hidden failure.

| Domain | Status | Evidence / boundary |
|---|---|---|
| Correctness | PASS | Full tests, official RFU integration, evidence index |
| Determinism | PASS | Restart/chunk/order properties and representative VDJdb reproduction |
| Scalability | PASS | Full Wells validation and versioned 10k/100k/500k synthetic baseline |
| Failure recovery | PASS | Injected R, output, checksum, cache, and interrupted-run failures |
| API stability | PASS | Machine-readable 0.4 snapshot with stable-removal guard |
| Documentation | PASS | Navigable, warning-clean local Sphinx API build and synthetic tutorial |
| Portability | PARTIAL | Cross-OS Python-only workflow is configured but cannot be observed until reviewed and pushed |
| CI | PARTIAL | Linux Python 3.10/3.11/3.12 and OS smoke jobs configured; current unpushed tree has no remote run |
| Installation | PASS | Isolated wheel and sdist import, CLI, fixture, and tutorial checks |
| Provenance | PASS | Hashed external evidence manifests and portable repository index |
| Privacy | PASS | Tracked files and distributions contain no public/private runtime datasets or external databases/assets |
| Dependency hygiene | PASS | Matplotlib/MuData remain optional; core clean install and import pass |
| TCR scientific validation | PASS | Official parity, Wells scale, independent GSE190905, preregistered GSE157007 |
| VDJdb external validation | PASS | Pinned 2026-06-03 reference; 24 analyses plus deterministic representative checks |
| BCR status | PARTIAL | Real-data preprocessing/features pass; frozen BCR reference correctly remains NO-GO |

## Release blockers

- Maintainer review of the complete uncommitted diff.
- A successful remote CI run on the exact reviewed commit, especially macOS and
  Windows smoke jobs, before promoting the candidate.
- Final release authorization and version promotion from `0.4.0rc1` to `0.4.0`.

## Strong recommendations, not blockers

- Publish the already buildable API documentation when release hosting is
  selected.
- Repeat the performance baseline on designated release hardware.
- Retest against the official external RFU checkout after any backend change.

BCR functional-unit construction, manuscript work, DOI creation, and additional
public cohorts are not required for the TCR software release.

