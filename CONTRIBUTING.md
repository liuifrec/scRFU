# Contributing to scRFU

Please open an issue before a large change. Bug reports should include the
scRFU version, Python version, operating system, a minimal synthetic example,
and `scrfu doctor` output. Do not attach private receptor data, credentials, RFU
assets, VDJdb, or large public datasets.

Development requires Python 3.10 or newer:

```bash
python -m pip install -e ".[dev,docs]"
python -m pytest -q
ruff check .
ruff format --check .
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

Contributions are submitted under the repository's MIT license. Follow the
project code of conduct and avoid including identifying or confidential data in
issues, logs, tests, or patches.

