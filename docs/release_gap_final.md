# Final public-release gap

Updated 2026-08-25. This is a release-readiness assessment, not approval to tag,
publish, archive, or upload. Existing Git tags are `v0.2.0` and `v0.3.0`; the
source is now prepared as **0.4.0**. This source version has not been tagged,
published, archived, or released.

## Current state

| Release item | State | Evidence | Remaining action |
|---|---|---|---|
| Package metadata | Complete | `pyproject.toml` has Yu-Chen Liu, `liu_y@rerf.or.jp`, and the canonical homepage/repository/issues URLs | Recheck built metadata in the release candidate. |
| Version reconciliation | Complete for preparation | source, citation metadata and tests use 0.4.0; historical tags remain untouched | Verify the exact eventual release commit before tagging. |
| CI | Complete for repository checks | `.github/workflows/ci.yml` tests Python 3.10–3.12, Ruff, build, wheel install, import, CLI, and a small installed-package smoke test | Require green CI on the exact release commit. |
| Public installation test | Strong partial | CI installs the locally built wheel into a clean environment | Test the final artifact from the intended public distribution channel after upload to a staging service or from the exact downloadable release artifact. |
| API documentation | Complete locally | README, warning-clean Sphinx/autodoc sources, the API snapshot, focused method docs, and examples exist | Hosted deployment remains optional until separately authorized. |
| Tutorial | Complete | `docs/tutorial.md` and `examples/tutorial_end_to_end.py` cover canonical receptors, mock/official backend distinction, metrics, longitudinal utilities, synthetic VDJdb and experimental BCR preprocessing | Keep the clean-wheel smoke green. |
| Small reproducible public test dataset | Complete | three frozen fully synthetic TSV fixtures with documented hashes are packaged | Preserve the hashes or version fixture changes explicitly. |
| `CITATION.cff` | Complete for preparation | authorship, repository, title, license and 0.4.0 are present without DOI/ORCID invention | Add a DOI only after one actually exists. |
| Changelog | Complete for preparation | `CHANGELOG.md` has a substantive 0.4.0 entry without inventing a release date | Preserve it on the eventual release commit. |
| License | Complete | `LICENSE` is present and package metadata identifies MIT | Include and inspect it in both sdist and wheel. |
| Zenodo-ready archive | Blocked | no approved versioned release or DOI | Exclude private/public source data, RFU assets, VDJdb, caches, and runtime results; archive only after the exact release artifact passes inspection. |
| Release notes | Complete for preparation | `docs/release_notes_0.4.0.md` separates package capabilities, external requirements, and limitations | Review with the final artifact hashes. |

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
2. Confirm source/package/changelog/release-note version is consistently 0.4.0.
3. Run the full test, formatting, build, wheel-install, integration, and artifact
   privacy inspection on the exact candidate.
4. Review the diff and obtain release approval.
5. Only then create a tag, public release, and archive/DOI.

None of these steps has been authorized or executed in this validation pass.
