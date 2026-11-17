# Versioning and release policy

The current software integration candidate is `0.4.1rc1` (unreleased). See
[current status](current_status.md) and the
[manual release checklist](release_notes_0.4.1rc1.md). Runtime, CITATION.cff,
wheel and sdist metadata must agree. Historical tags remain unchanged. A source
version change does not create a tag, package upload, hosted documentation,
GitHub Release, archive or DOI.

## Historical 0.4.0 preparation policy

The following text records the earlier release-preparation checkpoint; it is
not a statement of current publication status.

The release-preparation source reports version `0.4.0`. Historical Git tags
`v0.2.0` and `v0.3.0` remain unchanged. Changing the source version does not by
itself create a tag, package upload, GitHub Release, archive, or DOI.

Before an approved 0.4.0 release, confirm that runtime, wheel, sdist, changelog,
release notes, and the future tag all report 0.4.0. The prepared source must not
be published or tagged until clean-install, integration, artifact, and
maintainer-review gates pass.

For every future release:

1. the runtime version and wheel/sdist metadata must match;
2. a Git tag `vX.Y.Z` must point to exactly that source version;
3. the changelog must name the same version and compatibility changes;
4. the GitHub Release and DOI archive, when created, must reference the same tag
   and immutable artifacts.

Before 1.0, public stable APIs follow semantic-version intent: incompatible
changes require a documented migration and appropriate minor-version increase;
compatibility aliases receive `FutureWarning` with a replacement and removal
version before removal. Experimental APIs can gain fields while their scientific
defaults are being validated, but existing fields are not silently
reinterpreted.
