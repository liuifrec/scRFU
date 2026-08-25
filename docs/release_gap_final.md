# Final public-release gap

Updated 2026-08-25. This is a release-readiness assessment, not approval to tag,
publish, archive, or upload. Existing Git tags are `v0.2.0` and `v0.3.0`; the
source is now the unreleased development candidate **0.4.0rc1**. The next final
public version remains **0.4.0** after release review.

## Current state

| Release item | State | Evidence | Remaining action |
|---|---|---|---|
| Package metadata | Complete | `pyproject.toml` has Yu-Chen Liu, `liu_y@rerf.or.jp`, and the canonical homepage/repository/issues URLs | Recheck built metadata in the release candidate. |
| Version reconciliation | Strong partial | source, citation metadata and tests use unreleased 0.4.0rc1; historical tags remain untouched | Change to final 0.4.0 only on explicit release approval. |
| CI | Complete for repository checks | `.github/workflows/ci.yml` tests Python 3.10–3.12, Ruff, build, wheel install, import, CLI, and a small installed-package smoke test | Require green CI on the exact release commit. |
| Public installation test | Strong partial | CI installs the locally built wheel into a clean environment | Test the final artifact from the intended public distribution channel after upload to a staging service or from the exact downloadable release artifact. |
| API documentation | Strong partial | README, `docs/api_freeze_v1.md`, focused method docs, and examples exist | Build and review one navigable public API site or explicitly choose README/MkDocs/Sphinx delivery; verify all public symbols and versioned examples. |
| Tutorial | Complete | `docs/tutorial.md` and `examples/tutorial_end_to_end.py` cover canonical receptors, mock/official backend distinction, metrics, longitudinal utilities, synthetic VDJdb and experimental BCR preprocessing | Keep the clean-wheel smoke green. |
| Small reproducible public test dataset | Complete | three frozen fully synthetic TSV fixtures with documented hashes are packaged | Preserve the hashes or version fixture changes explicitly. |
| `CITATION.cff` | Complete for candidate | authorship, repository, title, license and 0.4.0rc1 are present without DOI/ORCID invention | Change only the version at final approval. |
| Changelog | Strong partial | `CHANGELOG.md` has a substantive Unreleased section | Convert Unreleased to 0.4.0 with the actual release date only at release approval. |
| License | Complete | `LICENSE` is present and package metadata identifies MIT | Include and inspect it in both sdist and wheel. |
| Zenodo-ready archive | Blocked | no approved versioned release or DOI | Exclude private/public source data, RFU assets, VDJdb, caches, and runtime results; archive only after the exact release artifact passes inspection. |
| Release notes | Weak partial | changelog can seed release notes | Draft concise 0.4.0 notes that separate demonstrated evidence from available APIs and name the external RFU requirement. |

## Scientific evidence relevant to the software release

- Public GSE190905 supplies the repeated-measures software validation; no
  private cohort is required for the independent tool release.
- Release-pinned VDJdb 2026-06-03 matching and null analyses are complete. The
  external database remains excluded from artifacts.
- The full receptor-only Wells assignment and downstream stress pass are
  complete; runtime outputs remain external.
- Final source-data deposits must contain only disclosure-reviewed derived
  public tables. No private cohort inputs or identifiers, VDJdb records, RFU
  assets, H5AD, or source-cohort databases may be bundled.
- Figure/source-table checksums and the final manuscript claim audit must be
  updated after either blocked analysis is executed.

## Minimal release sequence after scientific approval

1. Complete public API/docstring review and resolve final audit blockers.
2. Set source/package/changelog/release-note version to 0.4.0 only after
   approval.
3. Run the full test, formatting, build, wheel-install, integration, and artifact
   privacy inspection on the exact candidate.
4. Review the diff and obtain release approval.
5. Only then create a tag, public release, and archive/DOI.

None of these steps has been authorized or executed in this validation pass.
