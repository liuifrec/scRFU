# Final public-release gap

Audit date: 2026-08-24. This is a release-readiness assessment, not approval to
tag, publish, archive, or upload. Existing Git tags are `v0.2.0` and `v0.3.0`,
while the package remains 0.1.0. The next coherent public version should
therefore be **0.4.0**, assigned only after the manuscript-scope validation and
release checklist are approved.

## Current state

| Release item | State | Evidence | Remaining action |
|---|---|---|---|
| Package metadata | Complete | `pyproject.toml` has Yu-Chen Liu, `liu_y@rerf.or.jp`, and the canonical homepage/repository/issues URLs | Recheck built metadata in the release candidate. |
| Version reconciliation | Blocked | source version is 0.1.0; historical tags are `v0.2.0` and `v0.3.0` | Set package/changelog/release notes coherently to 0.4.0 only after approval; do not rewrite historical tags. |
| CI | Complete for repository checks | `.github/workflows/ci.yml` tests Python 3.10–3.12, Ruff, build, wheel install, import, CLI, and a small installed-package smoke test | Require green CI on the exact release commit. |
| Public installation test | Strong partial | CI installs the locally built wheel into a clean environment | Test the final artifact from the intended public distribution channel after upload to a staging service or from the exact downloadable release artifact. |
| API documentation | Strong partial | README, `docs/api_freeze_v1.md`, focused method docs, and examples exist | Build and review one navigable public API site or explicitly choose README/MkDocs/Sphinx delivery; verify all public symbols and versioned examples. |
| Tutorial | Weak partial | focused examples and `examples/README.md` exist | Write one end-to-end public TCR tutorial using a redistributable small input and external RFU setup, then execute it from a clean wheel. |
| Small reproducible public test dataset | Weak partial | synthetic fixtures support offline tests; public cohort manifests support acquisition | Choose and license-review a tiny redistributable receptor fixture, or clearly designate a generated synthetic dataset as the official public tutorial dataset with deterministic generator and expected hashes. Do not commit source H5AD, RFU assets, or restricted data. |
| `CITATION.cff` | Missing | no file is tracked | Add project authorship, repository URL, title, license, and version after the release version is frozen; do not invent a DOI or ORCID. |
| Changelog | Strong partial | `CHANGELOG.md` has a substantive Unreleased section | Convert Unreleased to 0.4.0 with the actual release date only at release approval. |
| License | Complete | `LICENSE` is present and package metadata identifies MIT | Include and inspect it in both sdist and wheel. |
| Zenodo-ready archive | Blocked | no approved versioned release or DOI | Exclude private/public source data, RFU assets, VDJdb, caches, and runtime results; archive only after the exact release artifact passes inspection. |
| Release notes | Weak partial | changelog can seed release notes | Draft concise 0.4.0 notes that separate demonstrated evidence from available APIs and name the external RFU requirement. |

## Scientific and manuscript prerequisites

- Governed longitudinal evidence is the minimum scientific blocker if the
  manuscript retains its longitudinal representation claim or Figure 2.
- A release-pinned VDJdb run is required only if the antigen-coherence Results
  claim and primary Figure 4 are retained. Otherwise use the prespecified
  robustness/transfer fallback and remove the antigen claim.
- Final source-data deposits must contain only disclosure-reviewed derived
  public tables. No private cohort inputs or identifiers, VDJdb records, RFU
  assets, H5AD, or source-cohort databases may be bundled.
- Figure/source-table checksums and the final manuscript claim audit must be
  updated after either blocked analysis is executed.

## Minimal release sequence after scientific approval

1. Decide whether longitudinal and VDJdb claims are included or removed.
2. Complete the tutorial, official small test fixture decision, and public API
   review.
3. Add `CITATION.cff` without speculative identifiers.
4. Set source/package/changelog/release-note version to 0.4.0.
5. Run the full test, formatting, build, wheel-install, integration, and artifact
   privacy inspection on the exact candidate.
6. Review the diff and obtain release approval.
7. Only then create a tag, public release, and archive/DOI.

None of these steps has been authorized or executed in this validation pass.
