# Versioning and release policy

The source package currently reports 0.1.0. Historical Git tags and package
metadata have not been demonstrated to be mutually consistent, so no existing
tag is treated as a methods-paper release solely by its name. This Month 1 work
does not change, move, or create tags and does not change the package version.

Before the next release, inspect each historical tag, confirm whether it is
annotated or lightweight, compare its package version, and document whether a
GitHub Release exists. The next version must not reuse a historical tag. Given
existing `v0.2.0` and `v0.3.0` names, 0.4.0 is the default recommendation unless
the audit establishes a different coherent policy.

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
