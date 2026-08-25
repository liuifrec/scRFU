# Release-candidate audit

Audit date: 2026-08-25. Candidate version: 0.4.0rc1. This document does not
authorize a commit, tag, upload, archive or publication.

## Release blockers

| Item | State | Required action |
|---|---|---|
| Exact candidate validation | Complete for current uncommitted tree | 324 tests collected: 320 passed and 4 optional integrations skipped; Ruff/build/official integration and isolated-wheel tutorial passed. Rerun on the exact approved release commit. |
| Artifact privacy/content audit | Complete for current artifacts | Wheel and sdist contain intended source/docs/examples/tests and three synthetic fixtures; no runtime paths, private data, caches, RFU/VDJdb assets or credentials were found. Repeat after any edit. |
| Release approval | Blocked by policy | A maintainer must review the uncommitted diff and explicitly approve changing 0.4.0rc1 to 0.4.0 and releasing. |

Private cohort analysis, manuscript formatting, a DOI and a BCR functional-unit
reference are not software release blockers.

## Strong recommendations

| Item | Current evidence | Recommendation |
|---|---|---|
| API reference presentation | Exported APIs have annotations/docstrings and focused method docs | Choose and review a navigable Sphinx/API site before or immediately after 0.4.0. |
| Manual external integration | Official RFU parity/integration is locally reproducible | Keep external integration documented and separate from core CI because RFU assets and R are not bundled. |
| Performance tracking | 10k/100k/500k synthetic baseline and advisory comparator exist | Re-run on release hardware; investigate warnings, never use fragile CI cutoffs. |
| BCR experimental surface | Canonicalization, pairing, state features and real GSE266519 adapter QC pass | Keep experimental labeling and NO-GO reference gate until donor/cohort holdout criteria pass. |

## Optional improvements

- Add hosted API documentation and example notebooks without making notebook
  execution part of core CI.
- Add another public repeated-measures TCR cohort if a future biological claim
  requires temporal depth beyond GSE190905.
- Add a BCR reference prototype only after the documented gate passes.
- Create a DOI/archive after, not before, an approved final release.

## Candidate contents

- Python 3.10/3.11/3.12 CI, package build, wheel import/CLI smoke and synthetic
  tutorial smoke.
- Lazy plotting and optional MuData extras; core import has no Matplotlib
  import.
- Synthetic packaged fixtures only; no RFU asset or real VDJdb row.
- Full Wells, public-cohort and VDJdb runtime evidence stays outside Git.
- `CITATION.cff`, MIT license, changelog, version policy and release audit.
