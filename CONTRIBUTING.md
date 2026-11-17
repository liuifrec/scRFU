# Contributing to scRFU

Please open an issue before a large change. Bug reports should include the
scRFU version, Python version, operating system, a minimal synthetic example,
and `scrfu doctor` output. Do not attach private receptor data, credentials, RFU
assets, VDJdb, or large public datasets.

Development requires Python 3.10 or newer:

```bash
python -m pip install ".[dev,docs,all]" openpyxl pyarrow
python tools/run_installed_tests.py --python-only
ruff check src tests examples tools docs/conf.py
ruff format --check src tests examples tools docs/conf.py
sphinx-build -W --keep-going -b html docs docs/_build/html
python -m build
```

Core tests must remain independent of R, RFU assets, VDJdb, and network access.
Integration changes should use a user-supplied external RFU checkout and record
artifact hashes without redistributing it. New tests and tutorials must use
small synthetic, fully redistributable fixtures.

Canonical TCR RFU definitions and frozen references require explicit scientific
review. BCR APIs are experimental preprocessing/feature utilities only; a BCR
functional-unit API must not be introduced unless the documented feasibility
gate is independently satisfied.

Contributions are submitted under the repository's MIT license. Avoid including identifying or confidential data in
issues, logs, tests, or patches.


Use `python tools/check_distributions.py` after building to verify separate clean
wheel/sdist installs outside the checkout. CI and local lint use the same software
paths above. The complete software suite includes five external-R checks in
`tests/test_integration_rfu_repo.py`; set `RFU_DIR` and, for the extension-only
check, `SCRFU_MAP_AWARE_RFU_DIR`. Without these, report skipped backend tests
explicitly. Standard CI runs Python tests separately from these external-R tests.
Manuscript validation belongs to its separate scientific workstream.
