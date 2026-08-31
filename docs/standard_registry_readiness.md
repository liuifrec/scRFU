# Standard-registry readiness

Audit date: 2026-08-31. No package was published.

- `python -m build` produces a universal wheel and source distribution.
- `twine check` passes for both current artifacts.
- PyPI returned no matching distribution for the exact project name `scrfu` at
  audit time. This is not a reservation and can change before publication.
- Clean wheel and sdist installs pass from outside the source tree, including
  the optional Scirpy/MuData tutorial.
- The package does not contain upstream RFU assets, VDJdb, or validation data.

For a future credential-free PyPI release, the maintainer should create a
protected GitHub environment named `pypi`, configure a PyPI trusted-publisher
record for `liuifrec/scRFU` and the reviewed release workflow/environment, and
use a tag-triggered workflow with `id-token: write` plus
`pypa/gh-action-pypi-publish`. The workflow must build and validate artifacts
before entering the protected publish environment. No API token should be
stored in the repository.

Conda-forge is feasible after a public source/PyPI release supplies a stable
URL and hash. A feedstock recipe should depend only on scRFU's Python package
dependencies and must not download or redistribute the external RFU checkout.
Bioconda is not required merely to satisfy the standard-registry criterion and
should be evaluated separately if a workflow-specific need appears.

