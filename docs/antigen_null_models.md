# Antigen-Coherence Null Models and Baselines

`scrfu.tl.rfu_antigen_permutation_test()` benchmarks an observed coherence
statistic against RFU-label permutations across eligible unique sequences. A
permutation changes labels only; it preserves the observed RFU size multiset,
the matched sequences, and their antigen annotations. The default statistic is
the within-RFU same-antigen pair fraction.

The unstratified size-preserving permutation is the primary null. Optional
strata are CDR3 amino-acid length, TRBV (with allele suffix removed), or both.
Within-stratum shuffling preserves RFU sizes globally and controls only the
listed sequence properties. It is not a causal model and does not establish
functional equivalence.

For observed statistic \(T\) and \(B\) permutations, the empirical upper-tail
probability is

\[
  p = \frac{1 + \sum_{b=1}^{B} 1[T_b \ge T]}{1+B}.
\]

The result also records the null mean, population standard deviation, a z score
when variance is nonzero, seed, stratification, eligible/matched sequence and
RFU counts, and the null values. A zero-variance null has an undefined z score.
The input is never mutated.

`compare_antigen_groupings()` provides deliberately simple descriptive
baselines: RFU, TRBV, exact CDR3 length, TRBV plus length, and one deterministic
size-matched random partition. It returns tidy purity, entropy, and pair
fraction rows. These baselines are not exhaustive clustering benchmarks, and
repeated exact CDR3 sequences are not presented as an independently discovered
grouping.

