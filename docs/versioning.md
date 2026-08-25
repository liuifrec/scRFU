# Versioning and release policy

The source package now reports the unreleased development version `0.4.0rc1`.
Historical Git tags `v0.2.0` and `v0.3.0` already exist, so the next coherent
public release line is 0.4.0. No historical tag is moved or reused, and this
development version is not itself a tag or published release.

Before the final 0.4.0 release, confirm that runtime, wheel, sdist, changelog,
release notes, and the future tag all report 0.4.0. The release candidate must
not be published or tagged until clean-install and integration gates pass.

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
