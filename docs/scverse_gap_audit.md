# scverse ecosystem gap audit

Audit date: 2026-09-01. The source of truth is the current scverse ecosystem
package registry checklist and schema, not historical project planning.

| Requirement | Status | Evidence or action |
|---|---|---|
| Public OSI-licensed code | PASS | Public GitHub repository; MIT `LICENSE` |
| Versioned releases | PASS | Existing version tags and release history |
| Standard-registry installation | BLOCKED | PyPI/conda publication has not occurred; source/wheel installation is validated |
| Essential automated tests | PASS | Unit, property, integration, serialization and clean-install tests |
| CI on push/PR | PASS | Exact commit `a6d48b60…998d` passed package, optional, docs, Python 3.10/3.11/3.12, Ubuntu, macOS, and Windows jobs in run `33362695818` |
| Public API documentation | PARTIAL | Sphinx sources and strict build pass; hosted service is not configured yet |
| Appropriate scverse data structures | PASS | Native AnnData/MuData AIRR modality, aligned Awkward chain records |
| Maintainer agreement to listing | BLOCKED | Requires explicit maintainer confirmation at submission time |
| Tutorials | PASS | General and scverse-native synthetic tutorials |
| Cookiecutter-scverse template | NOT REQUIRED | Recommended, not mandatory; current mature repository is not template-derived |

## External manual actions

1. Maintainer reviews and publishes a versioned package to a standard registry
   (PyPI is the prepared first target).
2. Maintainer provisions hosted documentation from `.readthedocs.yaml` and
   verifies that the public build exposes `objects.inv`.
3. Release CI is confirmed for the exact published commit and stable package
   and documentation URLs replace placeholders in the registry draft.
4. Maintainer explicitly agrees to the scverse code of conduct and ecosystem
   listing requirements.
5. Maintainer reviews `scverse_meta_draft.yaml` against the then-current schema
   and submits the registry pull request.

These are external administrative/release actions, not missing scientific or
local engineering functionality. No package, documentation service, or
registry PR was published or created by this work.
