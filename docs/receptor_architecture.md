# Receptor Architecture

## Audit and compatibility policy

Before this refactor, `wells.py` mixed generic encoded-dataframe HDF5 reading
with Wells cache policy; `extract.py` combined AIRR aliases with special
`uns["TCR_IR"]` handling; and the Wells example implemented VDJ slot selection
and result writing itself. The RFU backend was already mostly generic, but its
public entry point began with AnnData extraction.

Genuinely Wells-specific behavior is now limited to the `wells_tcr_ir` adapter
and legacy wrappers: locating `TCR_IR`, interpreting flattened VDJ 1/2 fields,
primary selection by duplicate count, consensus count, then slot order, and
recognizing the public-atlas layout.

Dataset-independent layers now cover selective H5AD reading, the canonical
schema and validator, adapters, the portable cache, table-level RFU execution,
metadata joins, and descriptive metrics.

```text
Wells TCR_IR ─┐
scirpy AIRR ──┼─> adapter -> canonical receptors + cell metadata
AIRR table ───┘                    |
                                  v
                         portable receptor cache
                                  |
                                  v
                     call_rfu_table (dataset-independent)
                                  |
                     per-sequence + per-row + mapping
```

Public Wells helpers remain compatibility APIs. New work should use
`prepare-receptors`, `read_receptor_cache()`, and `call_rfu_table()`. Legacy
cache migration is explicit and non-overwriting by default. Wells remains an
adapter, scalability benchmark, and case study—not the package data model.
