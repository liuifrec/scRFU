# RFU Phenotype Coupling

`scrfu.tl.rfu_phenotype_coupling()` returns one row per RFU and observed
phenotype level. It reports total and phenotype-specific abundance, proportion,
dominant phenotype/fraction, represented phenotype count, and entropy.

For RFU phenotype proportions `p_k`, entropy is `H=-sum(p_k log(p_k))`.
Normalized entropy is `H/log(K)`, where K is the number of phenotype levels in
the supplied dataset (zero when K=1). Phenotype specificity is exactly
`1-normalized_entropy`: one denotes concentration in one phenotype and zero an
even distribution across all supplied phenotypes. This is descriptive, not a
test of association.

With `sample_key`, overall recurrence is the number of independent samples
containing the RFU and overall prevalence divides by all samples. The long
table separately reports recurrence and prevalence within each phenotype.
Cell abundance is therefore never mislabeled as donor/sample recurrence.
Nearest/threshold assignment and cell/unique-sequence weighting are explicit;
rare RFUs are kept subject only to the stated `min_abundance`.
