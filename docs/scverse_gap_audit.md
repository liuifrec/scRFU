# scverse ecosystem gap audit

Audit date: 2026-08-31. The source of truth is the current scverse ecosystem
package registry checklist and schema, not historical project planning.

| Requirement | Status | Evidence or action |
|---|---|---|
| Public OSI-licensed code | PASS | Public GitHub repository; MIT `LICENSE` |
| Versioned releases | PASS | Existing version tags and release history |
| Standard-registry installation | BLOCKED | PyPI/conda publication has not occurred; source/wheel installation is validated |
| Essential automated tests | PASS | Unit, property, integration, serialization and clean-install tests |
| CI on push/PR | PASS | Linux Python matrix and cross-platform Python smoke workflows |
| Public API documentation | PARTIAL | Sphinx sources and strict build pass; hosted service is not configured yet |
| Appropriate scverse data structures | PASS | Native AnnData/MuData AIRR modality, aligned Awkward chain records |
| Maintainer agreement to listing | BLOCKED | Requires explicit maintainer confirmation at submission time |
| Tutorials | PASS | General and scverse-native synthetic tutorials |
| Cookiecutter-scverse template | NOT REQUIRED | Recommended, not mandatory; current mature repository is not template-derived |

Additional actionable readiness gaps are hosted documentation setup, a real
standard-registry release, stable documentation URLs for the registry draft,
and maintainer review of the eventual ecosystem pull request. No registry PR is
submitted by this work.

