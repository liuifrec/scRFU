# RFU Pseudobulk

`scrfu.tl.rfu_pseudobulk()` creates a sample-by-RFU matrix. `sample_key` is
mandatory: every row is one biological sample/donor/library named by the user.
It never pools all cells from a condition into one pseudo-sample.

`weighting="cell"` counts unique cell-RFU units within sample;
`unique_sequence` counts unique exact-CDR3-RFU units. `nearest` uses all
non-missing nearest RFU assignments; `threshold_pass` first requires
`pass_thr=True`. Minimum cell/sequence requirements are recorded per sample in
the QC table. Prevalence filtering may be an integer number of samples or a
fraction of samples. RFU and sample order are deterministic.
Categorical zero-count RFU levels are retained only with
`retain_zero_rfus=True`; that explicit request takes precedence over prevalence
filtering for those zero levels.

Normalization modes are raw `count`, row `proportion`,
`counts_per_1000 = 1000 count / sample total`, and CLR. For RFU count `x_ij`,
CLR is `log(x_ij + p) - mean_j(log(x_ij + p))`, where positive pseudocount `p`
defaults to 0.5. Thus an all-zero row transforms to zeros. This is explicit
zero handling, not imputation of biological observations.

The `RFUPseudobulkResult` contains `matrix`, consistent sample metadata,
RFU abundance/prevalence metadata, per-sample QC, and parameters. If requested
phenotype metadata map a sample to conflicting labels, construction raises
instead of choosing a label.
