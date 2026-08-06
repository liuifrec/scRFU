# Descriptive RFU Metrics and Weighting Semantics

`scrfu.tl.rfu_metrics()` reports descriptive summaries for existing RFU
assignments. It does not fit a model, test a hypothesis, or alter an RFU call.
Both phenotype grouping and weighting must be supplied explicitly:

```python
metrics = scrfu.tl.rfu_metrics(
    per_cell_results,
    groupby=["cell_type", "tissue"],
    weighting="cell",  # or "unique_sequence"
    donor_col="donor_id",
    sample_col="library_id",
)
```

One output row represents one RFU within one explicitly named phenotype group.
Rows without an RFU or CDR3 are excluded. The default identifiers are
`cell_id`, `cdr3aa`, `rfu_label`, and `pass_thr`; callers may name alternatives.

## Definitions

Within phenotype group *g* and RFU *r*:

- **RFU cell count** is the number of distinct cells assigned to *r*.
- **RFU cell abundance** is that cell count divided by the number of distinct
  RFU-assigned cells in *g*. A cell assigned to multiple RFUs can contribute
  once to each RFU, so RFU abundances need not sum to one in multi-chain input.
- **Unique-CDR3 richness** is the number of distinct exact amino-acid CDR3
  sequences assigned to *r*.
- **RFU sequence convergence ratio** is the unique-CDR3 richness of *r*
  divided by the number of distinct RFU-assigned CDR3 sequences in *g*. It is
  the fraction of the group's observed sequence repertoire converging on that
  RFU definition; it is descriptive and does not imply convergent evolution.
- **Multiplicity** is the number of distinct cell–CDR3 observations assigned to
  *r*, divided by its unique-CDR3 richness. A value of one means every observed
  sequence occurs in one cell; values above one indicate repeated clonotypes.
- **Clonotype entropy** is Shannon entropy in natural-log units over cell counts
  of exact CDR3 clonotypes within *r*. It is zero for a single observed
  clonotype and increases as cell abundance is distributed across clonotypes.
- **Dominant-clonotype fraction** is the largest exact-CDR3 cell count divided
  by all distinct cell–CDR3 observations in *r*.
- **RFU threshold-pass rate** is the fraction of non-missing threshold statuses
  that pass, using the requested cell or unique-sequence unit.
- **Donor/sample prevalence** is the number of non-missing donors or samples
  containing at least one observation of *r*, divided by the number of
  non-missing donors or samples represented in *g*. Counts and denominators
  are reported with the prevalence.

## Cell versus unique-sequence weighting

`weighting="cell"` makes `weighted_abundance` the RFU cell abundance and
computes threshold-pass rate from one record per cell. Expanded clonotypes have
more influence. This answers a cell-composition question.

`weighting="unique_sequence"` makes `weighted_abundance` equal the sequence
convergence ratio and computes threshold-pass rate from one record per exact
CDR3. Each observed sequence has equal influence regardless of clonal size.
This answers a repertoire-diversity question.

Both named abundance columns, richness, multiplicity, entropy, and dominant
fraction are always returned so a manuscript can state exactly which unit was
used. Phenotype groups are never inferred from clustering fields or labels.

## Refined convergence names and compatibility

`cell_abundance` is the distinct-cell count. `convergence_richness` is the
distinct exact-CDR3 count within the RFU and intentionally equals
`unique_cdr3_richness`; it is a semantic alias, not another formula.
`normalized_convergence` is RFU exact-CDR3 richness divided by group-wide
exact-CDR3 richness. `mean_sequence_multiplicity` is distinct cell-CDR3
observations divided by RFU exact-CDR3 richness. `dominant_sequence_fraction`
is the largest exact-CDR3 cell count divided by all distinct cell-CDR3
observations. `threshold_pass_rate` uses the explicitly requested weighting.

Historical columns remain for compatibility: `rfu_cell_count` equals
`cell_abundance`; `rfu_cell_abundance` is a within-group cell proportion;
`sequence_convergence_ratio` equals `normalized_convergence`; `multiplicity`
equals `mean_sequence_multiplicity`; `dominant_clonotype_fraction` equals
`dominant_sequence_fraction`; and `rfu_threshold_pass_rate` equals
`threshold_pass_rate`. The older “clonotype” name refers to exact CDR3 here and
must not be interpreted as a supplied true clonotype identifier.

`assignment_policy="nearest"` includes every non-missing nearest assignment;
`assignment_policy="threshold_pass"` restricts metrics to threshold-qualified
assignments. An explicit `chain` filters canonical result tables through their
chain column. These controls are dataset-independent.

## Why B-cell ALFU is separate future work

RFU is defined here for T-cell receptor sequence representations and the
validated upstream RFU atlas. B-cell receptors undergo somatic hypermutation,
class switching, lineage evolution, and heavy/light-chain pairing constraints
that are not represented by these TCR assignments. A B-cell ALFU model therefore
requires its own biological definition, training/reference data, validation,
and weighting decisions. Relabeling TCR RFUs or reusing these descriptive
metrics would not create a validated B-cell model, so ALFU remains separate
future work.
