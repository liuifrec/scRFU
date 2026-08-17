# Conventional Repertoire Metrics

`scrfu.tl.repertoire_metrics(data, groupby=..., weighting=...)` accepts a
canonical receptor table, RFU per-row table, or AnnData observation table. It
is descriptive: one output row is one explicit phenotype/sample group, not one
statistical replicate inferred from cells.

Let abundance counts over identities be `n_i`, `N = sum(n_i)`, and
`p_i = n_i/N`. Shannon entropy is `-sum(p_i log(p_i))`; Simpson diversity is
`1-sum(p_i^2)`; inverse Simpson diversity is `1/sum(p_i^2)`. For sorted
abundances `x_(i)` and `m` identities, Gini clonality is
`2 sum(i x_(i))/(m sum(x_i)) - (m+1)/m`, ranging from zero for equal abundance
to larger values for unequal abundance. Dominant fractions are
`max(n_i)/N` for the named CDR3 or true clonotype identity.

The output also reports receptor-bearing source rows, unique cells, exact-CDR3
richness, true clonotype richness when supplied, explicitly named
sequence-based fallback richness, mean/median CDR3 amino-acid length,
productive fraction (productive non-missing rows divided by all rows with a
non-missing productive status), V-gene richness, and top-V-gene fraction
(`max V abundance / total V abundance`). Mean and median lengths use the number
of amino-acid characters in each receptor-bearing row. Missing true clonotype
IDs are never silently filled:
`diversity_identity="cdr3aa_sequence"`, `clonotype_fallback=True`, and
`clonotype_richness` is missing.

Weighting must be explicit. `cell` counts unique cells per identity;
`unique_sequence` counts unique CDR3s per identity; `clonotype` counts true
clonotype units and requires no claim that cells are independent replicates.
Empty inputs return an empty table with the stable output schema. No p-values
or inferential differential-abundance results are calculated.
