# Software integration verification

Validation is performed on this reconstructed software tree, independently of
the prior 590-case local manuscript/development candidate. Exact local results
and remote checks are recorded for the integration PR head. See
[scope and excluded study tests](software_integration_scope.md).

Local Linux verification of the reconstructed integration tree:

| Check | Actual result |
| --- | --- |
| Installed wheel, Python 3.12.14 / pandas 3.1.1 | 517 passed, zero skips; 257 dependency warnings retained |
| Installed sdist, Python 3.10.20 / pandas 2.3.3 | 517 passed, zero skips; 32 dependency warnings retained |
| Separate clean core-only wheel/sdist environments | Installed-path assertions, CLI/doctor, fixtures, mock walkthrough, wrapper overrides and regulatory exports passed |
| Real upstream RFU engineering fixture, R 4.1.2 | 10/10 eligible rows match among 12 reordered inputs; maximum score error 0.0 at tolerance 1e-12; zero threshold/reconstruction mismatches |
| Real tutorial and serial/thread/process wrapper checks | Passed on Linux with external assets |
| Software lint/format, metadata/CFF, dependency consistency | Passed |
| Sphinx with warnings as errors; rendered local links | Passed; 57 HTML pages, no broken local links |

The reference used for real parity is upstream commit
`ff7a1ea6aca28444a7ff8ec7ef09745e7341f99d`. Core Python package jobs on
macOS/Windows do not execute this R backend. Their actual results and the Python
3.11 job are recorded in the integration PR checks, not inferred from Linux.

Commands:

```bash
python -m build
python -m twine check dist/*
python tools/check_distributions.py
python tools/run_installed_tests.py
ruff check src tests examples tools docs/conf.py
ruff format --check src tests examples tools docs/conf.py
sphinx-build -E -W --keep-going -b html docs docs/_build/html
```

For the complete 517-case suite configure `RFU_DIR` and
`SCRFU_MAP_AWARE_RFU_DIR`; the latter is an optional, separately supplied
extension, not part of official RFU. Standard CI runs the 512 Python cases with
`--ignore=tests/test_integration_rfu_repo.py`. This separation does not certify
actual R execution on macOS or Windows. Core smoke excludes the native scverse
tutorial test because the separate optional job installs and exercises it.

Installed wheel and sdist checks use separate fresh environments and an external
working directory, unset PYTHONPATH, assert the package's installed path, and
verify CLI/doctor, exports, fixtures, tutorial outputs, wrapper selection and
regulatory CLI/API behavior. Optional dependencies remain optional.

No biological cohort, QTL search, prediction model, expression test, sampling
analysis or notebook is run as part of this software integration.
