# scRFU 0.4.1rc1 — unreleased software candidate

The integrated source restores the complete intended public API, packages the
scRFU-owned R wrapper and synthetic walkthrough, resolves wrapper defaults from
installed resources through parallel execution, and handles regulatory missing
identifiers consistently under pandas 2 and 3. Explicit wrapper overrides remain
supported and doctor reports the actual selection. Upstream RFU implementation
and model files remain external. No API was removed or replaced by a stub.

Core-only installation does not require plotting, Scirpy, MuData, R or model
files. Antigen-set intersection preserves Python set semantics under pandas 3.
See the [integration scope](software_integration_scope.md), [verification](software_verification.md)
and [current status](current_status.md). This integration contains no manuscript
revision, biological recomputation or public-data redistribution.

## Later manual release steps

1. Review the integrated source, complete required CI and repository review, and
   verify an installed copy of the merged public main. Keep immutable installation
   pins tied to an already reachable repaired source commit.
2. Decide the final release version explicitly; synchronize runtime, wheel/sdist,
   CITATION.cff and changelog. Preserve historical tags.
3. Re-run build, distribution, metadata, software tests/lint and warning-clean docs
   checks for the exact release candidate. Review dependency/platform/backend
   limitations and external-resource terms.
4. Tag, upload a package, create a release/archive DOI, deploy documentation or
   submit a manuscript only after separate authorization. No automatic publishing
   or documentation-deployment workflow is enabled here.
