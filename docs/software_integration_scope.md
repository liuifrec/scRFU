# Public software integration scope

This software was reconstructed on public main. Its ancestry does not include
unpublished manuscript commits. The complete package implementation is retained;
no intended public symbol is removed or replaced by a stub. Publication-specific
genetic-to-single-cell RFU identity is not asserted.

| Included | Boundary |
| --- | --- |
| All runtime modules, packaged scRFU-owned R wrapper, synthetic fixtures | Upstream RFU code/models are external |
| Generic receptor/cellular/longitudinal/regulatory APIs and adapters | No bundled study inputs or results |
| Synthetic tutorials, generic CLI examples, bounded real-RFU parity fixture | Mock labels are explicit; external-data runs are opt-in |
| API/metadata/regression tests, wheel/sdist gates, software CI | No study endpoints or notebooks execute |
| User/API documentation, MIT license, citation and manual release notes | No publishing, release or documentation deployment workflow |

The local scientific workstream retains manuscript drafts, study notebooks,
source data, cohort outputs, figure assemblies and author-review bundles. The
full local file inclusion/exclusion manifest records source hashes and individual
decisions. This public [machine-readable scope](software_integration_scope.json)
records the test boundary without including those scientific artifacts.

The source suite had 590 cases. Thirteen study-specific test files contribute
73 cases; one submission-checklist test also depends on scientific evidence
records. Those 74 remain in the scientific workstream. The candidate retains
516 source software cases and adds one distribution-scope check, for 517 cases.
All regulatory normalization, reference-QC and preparation regressions remain,
as do synthetic evidence-tier, unresolved-direction, API export and optional
integration tests. No software failure was hidden by changing this boundary.

Five retained cases require configured external R/reference assets. Linux local
validation runs them with real upstream assets and the separately identified
map-aware extension. Standard CI explicitly excludes that file from its Python
suite; those five are NOT RUN remotely, not counted as passes. Platform package
jobs verify Python imports, resources and mock subprocess behavior, not real R.

The existing public historical plans remain at their public-base contents,
with archival context in the navigation; one broken historical-plan link is
identified as belonging to the separate workstream. Newer scientific edits are excluded.
API snapshots labelled 0.4/0.5 record development contracts; the actual package
version is now 0.4.1; 0.4.1rc1 identifies the preceding integration checkpoint. Frozen manuscript/notebook lint findings belong to the
separate workstream. Both local and CI software lint cover only `src`, `tests`,
`examples`, `tools` and `docs/conf.py`.
