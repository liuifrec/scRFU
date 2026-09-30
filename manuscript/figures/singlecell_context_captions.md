# Wells single-cell context captions

These supporting assets supplement the quantitative benchmark figures. The frozen
Wells atlas supplies all three data types. Local GSE190905 files contain receptor
assignments and metadata without the expression matrix or an embedding;
GSE280982 likewise has no available expression matrix. No new data were downloaded.

## A — RFUs in the existing transcriptomic atlas

**Poster caption:** The frozen Wells UMAP places qualified assignments and three
donor-supported RFUs in source-annotated transcriptomic space; spatial overlap
does not establish cell-state or antigen specificity.

**Manuscript legend:** Coordinates are copied without fitting, rotation or
rescaling from `obsm/X_umap` in the frozen Wells atlas. The same 40,516 cells
appear in every panel: a deterministic 40,000-cell atlas background plus every
qualified cell belonging to RFU3526, RFU2114 or RFU527 (562 selected-RFU cells;
46 already in the background). Source cell-type annotations are unchanged;
shortened labels are mapped to the full annotations in `A_annotation_key.tsv`.
All 233,913 threshold-qualified primary-TRB assignments among 610,429 atlas cells
contribute to the displayed assignment count; the point cloud is the documented
display subset. Purple marks qualified cells in the assignment panel or the
named RFU in each highlight panel; gray marks the remaining displayed cells.
The RFUs were selected using donor support alone, ranked by donors with at least
10 qualified cells, then donors with at least five, total cells, and literal RFU
label. RFU3526, RFU2114 and RFU527 have 204, 196 and 162 qualified cells,
respectively. Highlighted-cell oversampling means these panels should not be
used to estimate atlas cell-type proportions or enrichment. No annotation,
embedding, receptor assignment, or transcriptomic clustering was recomputed.

**Take-home:** Frozen transcriptomic coordinates provide a reusable context for
RFU assignments without establishing cell-state specificity.

## B — Observed amino-acid sequence families

**Poster caption:** Three RFUs summarize observed CDR3 amino-acid families;
logos weight each distinct sequence once and show only the modal sequence length.

**Manuscript legend:** Logos use the actual CDR3 amino-acid strings from the
frozen, threshold-qualified Wells primary-TRB assignments (score ≥0.6). Each
distinct amino-acid string receives equal weight, regardless of cell abundance
or V call. Within each RFU, the modal length among distinct strings is selected
(shorter length breaks a tie); sequences of other lengths are retained in the
source table and are not aligned or padded into the logo. Letter height is the
observed per-position amino-acid frequency, not information content. The heading
gives the number of strings in the logo and their length; the support line gives
all qualified cells, all distinct amino-acid strings and all supporting donors
for the family. These denominators differ deliberately. Full length distributions,
individual strings, V calls and cell provenance are supplied in the source
tables. Colors distinguish broad residue classes for readability. Logos describe
sequence families, not consensus antigens, antigen specificity or function.

**Take-home:** RFUs can be displayed through their observed constituent sequences;
the modal-length logo is a partial sequence-family summary.

## C — Exploratory matched RFU expression contrast

**Poster caption:** In Wells, RFU3526 versus matched other qualified RFUs yields
no q < 0.05 genes across 12 paired donors; only 5–26 cells per donor/group support
this exploratory comparison.

**Manuscript legend:** This is a new, bounded figure-specific exploratory
contrast derived from an existing frozen expression object, not a frozen
radiotherapy result. RFU3526-positive cells were matched without replacement to
cells assigned to other qualified RFUs within each exact donor × tissue ×
source cell type × library × V-call stratum. Within each stratum, equal numbers
were retained using SHA256-based cell ordering, independently of expression.
Donors with at least five matched cells per group were included: 12 donors,
122 cells per group, 5–26 cells per donor/group. Counts from `raw/X` were verified
to be nonnegative integers, summed within donor/group and converted to
log₂(CPM + 1), with library size defined by all raw genes. The horizontal axis
is the equal-donor mean paired difference, RFU3526 minus matched background;
it is not a fitted count-model log fold change. Two-sided paired t-tests treat
donors as replicates. Genes require CPM ≥1 in at least 12 of 24 pseudobulks,
positive variance of paired differences (numerical tolerance 1e-14), and a false
source `feature_is_filtered` flag. TCR genes matching `^TR[ABDG][VDJC]` are excluded
to avoid directly testing assignment-defining receptor features. Benjamini–Hochberg
adjustment covers all 6,650 tested genes. The vertical axis shows −log₁₀(raw P)
to display the distribution when adjusted P values are uniformly high; point
color encodes q < 0.05, and the figure explicitly reports zero significant genes.
The minimum adjusted P is 0.9754. The four labels are the lowest adjusted
P values, with raw P value, absolute effect size and gene ID resolving ties;
none is statistically significant at q < 0.05. No contrast or gene filter was
tuned after observing the results. Small pseudobulks, clone dependence, residual
state differences and observational RFU membership limit interpretation; this
is neither causal RFU biology nor pre/post-radiotherapy DE. No detected
significance is not evidence of expression equivalence.

**Take-home:** This small, donor-paired exploratory contrast does not support a
transcriptomic difference at q < 0.05.

## D — Supporting single-cell context block

**Poster caption:** A frozen Wells UMAP, observed RFU sequence families and an
explicitly exploratory expression comparison provide supporting single-cell context.

**Manuscript legend:** This assembly reuses the RFU3526 highlight from A, all
three modal-length logos from B, and the single expression contrast from C.
Coordinates, sequences, effect sizes and adjusted significance are identical
to the standalone panels. The lower-right whitespace is reserved for later
composition beside BioRender artwork; it contains no missing measurement.
Refer to A–C for sampling, sequence-length and small-pseudobulk limitations.
The block complements rather than replaces the longitudinal, matched-control,
external-transportability and regulatory quantitative figures.

**Take-home:** Transcriptomic context and observed sequence families are supporting
views; they do not supply antigen or causal mechanistic evidence.

## Methods references

Donor replication motivates pseudobulk summaries; see
[Squair et al., Nature Communications (2021)](https://www.nature.com/articles/s41467-021-25960-2).
The implemented test and correction are documented in
[SciPy paired t-test](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ttest_rel.html)
and [SciPy false-discovery control](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.false_discovery_control.html).
This simple paired log-CPM analysis is exploratory; it is not an edgeR/DESeq2 or
limma-voom model.
