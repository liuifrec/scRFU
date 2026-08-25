# Performance regression checks

`examples/synthetic_scaling_benchmark.py` measures deterministic pure-Python
validation, deduplication, reconstruction, repertoire metrics, pseudobulk,
overlap and phenotype coupling at declared sizes. It writes only a compact TSV
and manifest; generated receptor rows are never persisted.

`benchmarks/synthetic_python_baseline.tsv` records the 2026-08-25 Python 3.10.20
baseline for 10k, 100k and 500k rows on Linux. Peak RSS is process-cumulative,
so it is useful for detecting large shifts but is not attributed exclusively to
the named operation.

Compare a new external run with:

```bash
python examples/compare_benchmark_baseline.py \
  --baseline benchmarks/synthetic_python_baseline.tsv \
  --current /path/to/synthetic_scaling.tsv \
  --output /path/to/performance_comparison.tsv
```

Ratios above 1.5 are advisory warnings. The tool deliberately returns success
because wall time and RSS vary by operating system, allocator and hardware;
ordinary CI must not use fragile performance cutoffs. Output dimensions should
still be reviewed as correctness checks.
