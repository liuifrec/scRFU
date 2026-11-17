# Reference-coverage diagnostics

Per-row outputs retain nearest assignment and threshold qualification
separately. Explicit aliases are `rfu_id_nearest`, `rfu_label_nearest`,
`rfu_score`, and `rfu_pass_threshold`. `reference_coverage_status` is one of:

- `covered`: eligible, assigned, and threshold qualified;
- `below_threshold`: eligible with a nearest assignment below threshold;
- `ineligible_sequence`: present sequence that fails the frozen eligibility
  rule;
- `upstream_unassigned`: eligible query with no returned assignment;
- `missing_sequence`: no usable amino-acid CDR3.

`assignment_status` further names threshold-qualified and below-threshold
nearest assignments. Legacy `rfu_id`, `rfu_label`, and `pass_thr` remain
compatibility fields.

`scrfu.tl.reference_coverage()` reports input, eligibility, assignment and status
counts; assignment-score quantiles; threshold-pass fraction among eligible
rows; unique-sequence coverage; and unique-cell-weighted coverage, optionally by
sample/cohort/chain.

Threshold failure means low similarity relative to the frozen reference under
the selected threshold. It is not a calibrated confidence interval, a
probability of correct biological function, or a validated out-of-distribution
classifier.
