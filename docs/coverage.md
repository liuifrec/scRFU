# Test coverage

Coverage was measured on 2026-08-31 with the complete locally installed
optional scverse stack and branch tests:

```text
366 tests collected
362 passed, 4 optional integration skips
5,154 statements; 4,364 covered
79.92% combined line/branch coverage report
84.67% statement coverage
67.29% branch coverage
```

Modules that implement stable public API entries have 85.86% weighted statement
coverage. The new native storage module is 89% covered and its primary user
paths—multi-chain extraction, explicit summaries, H5AD/H5MU round trips,
subsetting, concatenation, X independence and current Scirpy interaction—are
tested directly.

The weakest stable modules are the completed-large-run validator (48.7%),
sequence matrix utilities (70.1%), repertoire conveniences (71.2%), and I/O
(72.3%). `wells.py` is 47.6% covered overall but its remaining gaps primarily
involve external/large-file branches; this sprint did not chase coverage by
rerunning full Wells. Coverage is therefore classified PARTIAL rather than
inflated to a release gate.

Run locally with:

```bash
python -m pytest --cov=scrfu --cov-report=term-missing
```
