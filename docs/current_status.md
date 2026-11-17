# Current software candidate

**0.4.1rc1 is unreleased.** This branch integrates the complete reviewed software
implementation onto public main, without importing unpublished manuscript history.
The {download}`source README <../README.md>` pins installation to the already
public repaired software commit `4bb6bc18cf7c762c2e15255aaebe5af7bc16df05`. No package-registry
release, software DOI or hosted documentation service is claimed.

Start with the [synthetic quickstart](tutorial.md), [first real RFU run](first_real_rfu.md),
[configuration](configuration.md), [single-cell guide](scverse_tutorial.md), and
[regulatory guide](regulatory_triangulation.md). See [integration scope](software_integration_scope.md)
and [release notes](release_notes_0.4.1rc1.md).

Core metadata supports Python 3.10+ and pandas >=2.3,<4. The software matrix tests
Python 3.10/3.11 with pandas 2.3.3 and Python 3.12 with pandas 3.1.1. Package jobs
build wheel and sdist and install each into separate core-only environments on
Linux, macOS and Windows, running outside the checkout. Optional scverse/plotting
integration tests are separate from the clean core-only checks. Python 3.13+ and
pandas 4 are not verified support combinations.

Actual integration results are recorded in [software verification](software_verification.md).
Python package smoke success on a platform does not establish actual R-backend
execution there. Upstream R code and reference files remain externally supplied.

Historical public plans in the archive retain their original context; their
readiness language does not establish release, manuscript approval or submission.
The scientific workstream, frozen study outputs and author declarations are not
part of this software integration. Reference-number equality does not establish
cross-dataset biological identity.
