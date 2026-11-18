# Current software status

**0.4.1** finalizes the public software integration at
`f320ede184f9b9741a28df54a2cc755248b46eb0` (previously 0.4.1rc1). No unpublished
manuscript history is included. The {download}`source README <../README.md>`
selects the `v0.4.1` tag published with the verified GitHub Release. The release
notes distinguish source/tag publication from package-registry publication,
a software archive DOI and hosted documentation; unverified services are not claimed.

Start with the [synthetic quickstart](tutorial.md), [first real RFU run](first_real_rfu.md),
[configuration](configuration.md), [single-cell guide](scverse_tutorial.md), and
[regulatory guide](regulatory_triangulation.md). See [integration scope](software_integration_scope.md)
and [release notes](release_notes_0.4.1.md).

Core metadata supports Python 3.10+ and pandas >=2.3,<4. The software matrix tests
Python 3.10/3.11 with pandas 2.3.3 and Python 3.12 with pandas 3.1.1. Package jobs
build wheel and sdist and install each into separate core-only environments on
Linux, macOS and Windows, running outside the checkout. Optional scverse/plotting
integration tests are separate from the clean core-only checks. Python 3.13+ and
pandas 4 are not verified support combinations.

Actual Linux upstream R execution was verified on the integration checkpoint;
release validation repeats the bounded engineering fixture, not biological analyses.
Actual integration results are recorded in [software verification](software_verification.md).
Python package smoke success on a platform does not establish actual R-backend
execution there. Upstream R code and reference files remain externally supplied.

Historical public plans in the archive retain their original context; their
readiness language does not establish release, manuscript approval or submission.
The scientific workstream, frozen study outputs and author declarations are not
part of this software integration. Reference-number equality does not establish
cross-dataset biological identity.
