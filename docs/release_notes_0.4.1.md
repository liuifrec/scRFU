# scRFU 0.4.1

This software release finalizes the reviewed public integration at f320ede
(0.4.1rc1). It restores the complete intended API, packages the scRFU-owned R
wrapper and synthetic walkthrough, keeps wrapper resources alive through serial
and parallel execution, and fixes regulatory missing identifiers under pandas 2
and 3. Explicit wrapper overrides remain supported; doctor reports the actual
selection. Optional plotting/scverse dependencies remain optional. The final
version change does not alter scientific defaults, thresholds or reference files.

## Scope, credit and support

The software author is Yu-Chen Liu, as documented in CITATION.cff and package
metadata. Cite the upstream RFU method by Hu, Pan, Reid, Tworoger and Li,
Nature Communications 15, 8171 (2024), doi:10.1038/s41467-024-52522-z.
Scirpy/scverse attribution is retained in the user guides. No manuscript author
list, unpublished manuscript citation or author identifier is inferred.

Python/package checks cover Linux, macOS and Windows. The tested dependency
matrix is Python 3.10/3.11 with pandas 2.3.3 and Python 3.12 with pandas 3.1.1.
Actual upstream R execution is validated on Linux with R 4.1.2 and externally
supplied implementation/model files. Mock-wrapper and platform package checks
are not real R validation; real R on macOS/Windows is NOT RUN. Python 3.13+
and pandas 4 are not verified support combinations.

The 517-case software suite includes five opt-in external-R cases; standard CI
runs the other 512. Seventy-four manuscript-specific/checklist cases remain in
the separate scientific workstream. See the [scope](software_integration_scope.md)
and [integration evidence](software_verification.md).

Unpublished manuscripts, study notebooks, biological/source data, cohort results,
local review bundles and upstream RFU model/code assets are excluded. Original
scientific execution versions remain historical facts. No assignments, cohort
endpoints, QTL searches, prediction, expression, sampling or notebooks were rerun.

## Manual publication procedure

Publication is authorized for this software release after the exact candidate
passes normal CI/review and installed-distribution gates. No automatic publisher
is enabled. Preserve historical tags and never replace published distributions.

1. Freeze the source commit, build wheel and sdist, inspect their members, run
   `twine check`, and verify separate clean installations outside the checkout.
   Run the installed software suites, Linux RFU engineering fixture, software
   lint/format and warning-clean documentation build. Record the source commit,
   results and artifact SHA256 hashes outside the source tree.
2. Merge the preparation PR only after required checks and review. Tag the exact
   verified source `v0.4.1`. Create the GitHub Release with the same verified
   wheel/sdist and checksum/provenance records; do not rebuild or replace them.
3. With authorized PyPI account access or a configured trusted publisher, upload
   those same files to the production registry. A missing public project does
   not establish ownership. Verify downloaded registry hashes and a fresh
   `scrfu==0.4.1` installation before advertising that registry command.
4. Check the existing software archive/integration for an automatically created
   record before any deposit. Without an established archive account/integration,
   leave this step pending. Record only returned version/concept DOIs and preserve
   their distinction; never invent an identifier or duplicate a deposit.
5. Publish the warning-clean built documentation to the configured destination
   when its project/account access is available. Do not invent a hosted URL or
   choose an unrelated hosting service. Keep the built documentation ready locally.
6. Verify public tag/release/download/citation/documentation links and install from
   the public release outside the checkout. Record any post-publication citation
   or documentation updates in subsequent commits; frozen artifacts stay unchanged.

The README selects source/tag installation. GitHub Release, package-registry,
archive DOI and hosted documentation are distinct outputs; final publication
records must report each actual outcome and any concrete remaining setup.
