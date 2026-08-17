# Dataset-independent Visualization

The plotting API consumes generic tables or structured analysis results:

- `rfu_metric_heatmap`: RFU rows by sample columns, deterministic top-N and
  optional row z-standardization;
- `rfu_overlap_heatmap`: square similarity/distance result with direction in
  its label;
- `rfu_convergence`: richness versus abundance/multiplicity, optional
  phenotype colors and deterministic top annotations;
- `rfu_phenotype_heatmap`: long coupling output pivoted to RFU by phenotype;
- `repertoire_metric_comparison`: descriptive sample-matched scatter of RFU
  and conventional metrics.

Every function accepts `ax=None`, returns a matplotlib `Axes`, never calls
`show`, rejects empty input, and leaves titles/labels customizable. Matplotlib
is loaded only on plotting use and is installed with `scrfu[plotting]`; seaborn
is not required.

The antigen-evidence layer adds `rfu_antigen_heatmap` (RFU by antigen),
`rfu_antigen_coherence` (matched richness versus purity/entropy),
`antigen_permutation_distribution` (observed statistic against the explicit
null), and `rfu_antigen_bubble` (matched abundance and within-RFU proportion).
Top-N selection breaks abundance ties lexically and is deterministic. These
plots visualize external evidence; only the permutation view carries an
explicit benchmark, and none labels an RFU as antigen-specific.

`scrfu.tl.rfu_sequence_matrix()` provides an optional logo foundation without
a plotting dependency. It returns the exact position-frequency matrix for
unique-sequence, cell, or true-clonotype weighting. Left, right, center, and
conserved-end layouts are visualization conventions, not formal alignment.
Conserved-end mode anchors an initial C and terminal F/W when both occur. Logo
motifs must not be interpreted as functional evidence without validation.
Logomaker integration remains deferred: the reproducible frequency matrix is
the stable foundation, avoiding a new mandatory plotting dependency.
