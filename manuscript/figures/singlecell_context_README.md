# Single-cell context reusable figure pack v1

All data come from the existing Wells public atlas and its completed primary-TRB
assignments. Existing GSE190905/GSE280982 local files lack expression matrices;
the alternative RFU-positive comparison is therefore used for the volcano.
Previous RP1-14, cross-application and APBC asset packs remain unchanged.

## Contents and placement

| Figure | Poster native size (inches) | Recommended use |
|---|---:|---|
| A_umap_context | 24 × 15 | Poster context; manuscript supplement |
| B_rfu_logos | 18 × 12 | Poster sequence-family context; manuscript supplement |
| C_exploratory_volcano | 13 × 11 | Supplement or optional poster backup |
| D_context_block | 26 × 18 | Optional supporting poster block with reserved whitespace |

Each is in `figures/poster/` and `figures/manuscript/`, in PDF, SVG and
300-dpi PNG. Poster text is at least 22 pt at native placement; reducing width
also reduces text. Manuscript labels are at least 8 pt. Dense scatter points
are rasterized at 300 dpi inside PDF/SVG, with vector text, axes and sequence
logos. Variable-height logo glyphs encode frequency and are data marks, not
fixed-size labels. White background, DejaVu Sans, CD4 blue/CD8 orange source
annotation families, and purple RFU highlights match the existing asset style.
`preview/contact_sheet.png` shows all four poster exports. Every figure has a
one-line caption, fuller legend and take-home in `manifests/captions.md`.

## Reproduction

From `/home/liuyuchen/Github/scRFU`, use the existing environment:

```bash
PY=/home/liuyuchen/data/scRFU_regulatory_20260909/environment/bin/python
"$PY" -m manuscript.scripts.singlecell_context_figures --prepare-only --out /tmp/scrfu-context-prepared-v1
"$PY" -m manuscript.scripts.singlecell_context_figures --prepared /tmp/scrfu-context-prepared-v1 --out /home/liuyuchen/data/scRFU_radiation_methods_20260916/results/singlecell_context_reusable_figures_v1
"$PY" -m manuscript.scripts.singlecell_context_figures --verify-only
```

Preparation requires an empty directory and performs the one authorized bounded
DE contrast once. **For figure revisions, reuse the completed pack itself with
`--prepared <completed-pack> --out <new-directory>`; do not prepare again.**
Rendering reads frozen plotting tables, verifies selected source rows, and
never repeats statistical fitting. It rejects nonempty destinations. Verification
does not rewrite any file. No assignment pipeline, embedding fit, prediction
model, QTL search or longitudinal endpoint is called. The complete source object
is never loaded into memory: only metadata/coordinates and 244 selected raw
expression rows are read. SHA256 input verification streams the complete file.

## Source tables

All tables are TSV. Full column lists are in `manifests/table_columns.json`.
Counts are integers, coordinates retain the exact float32 source values through
float64 TSV serialization, and statistical values use 17 significant digits.

| Table | Rows / interpretation |
|---|---|
| A_umap_cells | Displayed cells, exact source row and original UMAP coordinates; donor alias, source cell type, tissue, assignment status, literal RFU and subset flags |
| A_annotation_assignment_counts | Entire atlas cell counts per original cell type and qualified / below-threshold / no-primary-TRB status |
| A_annotation_key | Full source label, shortened figure label and display color; no new cell-state calls |
| rfu_selection_ranking | All qualified RFUs ranked using donor/cell support only; records ≥5/≥10-cell donor counts and largest donor contribution |
| B_assigned_cells | Every qualified cell in the three selected RFUs, its actual amino-acid string, V call, donor/tissue/source type and frozen assignment score |
| B_sequences | One distinct amino-acid string per RFU; full cell/donor support, observed V calls, length and modal-length inclusion flag |
| B_logo_frequencies | RFU, position, amino acid, number of unique strings, denominator, frequency; 20 residues including zeros at every position |
| B_rfu_support | Family cell/donor/unique-AA counts, modal length, logo sequence/cell support, maximum donor cell count |
| C_matching_strata | Exact matching keys, case/background availability, retained count per group and donor inclusion; includes strata without a background match |
| C_matched_cells | All initially matched cells with source row, literal RFU, exact matching metadata and donor-level DE inclusion; raw counts are blank for excluded cells, never zero-filled |
| C_pseudobulk_samples | Donor × group pseudobulk identifier, selected cell count and total raw counts |
| C_pseudobulk_counts | Raw gene ID and count sum for every donor/group; all raw genes are retained |
| C_pseudobulk_log2cpm | Same genes and columns after CPM then log₂(CPM+1) |
| C_differential_expression | All raw genes with exact contrast label, effect, donor count, filter status, raw P, BH-adjusted P, −log₁₀(q), significance and label flags; volcano y is −log₁₀(raw P), point color/count report adjusted significance; untested P/q are blank |

Donor aliases use the same sorted-atlas mapping W01–W24 as the earlier Wells
summary. Exact source donor IDs are retained only in external matching tables.
Family logos count distinct AA strings; this is different from nucleotide
CDR3/V/J receptor identity used by the longitudinal benchmark figures.

## Fixed analysis and interpretation

RFU3526, RFU2114 and RFU527 maximize broad donor support under the declared
ranking. They were selected without inspecting expression differences or UMAP
positions. Modal-length subsets use one vote per distinct CDR3 amino-acid string;
no alignment, trimming or antigen interpretation is introduced.

The RFU3526 volcano is a **new, figure-specific exploratory donor-paired contrast**,
not an existing radiotherapy result. Within donor/tissue/source-cell-type/library/
V-call strata, deterministic case/background sampling without replacement gives
12 donors with ≥5 cells per group (122/group overall, 5–26 per donor/group).
Paired t-tests use log₂(CPM+1) donor pseudobulks; donor weights are equal. The
effect is a paired log-CPM difference, not a count-model fold-change coefficient.
Low-expression, source-filtered, zero-variance and TCR genes are excluded before
BH correction. **0/6,650 tested genes have q < 0.05 (minimum q = 0.9754).** The
vertical axis uses raw P to show the spread; adjusted significance is explicit
in color and the zero-hit count. Labels identify four top
ranked genes, not significant hits. Small cell numbers and residual state/clone
confounding preclude causal or RFU-specific claims. No detected significance
does not establish equivalence. Full details and method links are in the captions.

## Traceability and validation

`manifests/preparation.json` records exact input paths/hashes, specification,
software versions, preparation code hash and all frozen table hashes.
`manifests/figure_manifest.json` records source tables per figure, dimensions,
minimum fonts, generating command, code hashes and export hashes. Per-export
`*_value_checks.json` verifies actual scatter coordinates and logo frequencies
against the tables. `completion.json` hashes every delivered file except itself.
Validation rejoins each plotted cell to the frozen object and assignment table,
checks all selected sequences and donor support, re-reads only selected raw rows
to verify count sums, checks the stored log-CPM transform/effects/BH arithmetic,
and validates file existence, rendering, page count, PNG dimensions and DPI.
Rechecking counts is validation, not a repeated biological workflow or DE fit.

Only code, safe metadata and documentation belong in Git; source tables,
sequences, cell IDs, processed expression and generated exports remain external.
