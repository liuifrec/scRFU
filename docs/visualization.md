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

`scrfu.tl.rfu_sequence_matrix()` provides an optional logo foundation without
a plotting dependency. It returns the exact position-frequency matrix for
unique-sequence, cell, or true-clonotype weighting. Left, right, center, and
conserved-end layouts are visualization conventions, not formal alignment.
Conserved-end mode anchors an initial C and terminal F/W when both occur. Logo
motifs must not be interpreted as functional evidence without validation.
