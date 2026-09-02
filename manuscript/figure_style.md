# scRFU figure style system

## Design goal

Figures should read as a coherent methods paper: precise, restrained and
legible at journal-column scale. Decoration never competes with the result, and
RFU is visually identifiable without being made larger, brighter or otherwise
favored over comparators.

## Provisional canvas sizes

These dimensions are working defaults, not journal commitments.

| Use | Width | Typical height |
|---|---:|---:|
| Single-column panel/figure | 89 mm / 3.50 in | 60–150 mm |
| 1.5-column figure | 135 mm / 5.31 in | 90–180 mm |
| Double-column main figure | 183 mm / 7.20 in | 120–220 mm |
| PNG preview | same physical size | 300 dpi |
| Vector master | same physical size | PDF with embedded text/fonts |

Main prototypes use double-column width. Final dimensions should be adjusted
only after the target journal's production specification is confirmed.

## Typography

- Typeface: DejaVu Sans for prototypes; later substitute a journal-compatible
  sans-serif only as a global change.
- Panel letters: 10 pt, bold, upper-left, aligned to the panel bounding box.
- Panel titles: 8 pt, semibold where supported.
- Axis titles: 7.5 pt.
- Tick labels: 6.5–7 pt.
- Legends and annotations: 6.5–7 pt; never below 6 pt at final size.
- Figure-level titles are omitted from submission artwork; working PNGs may
  include a short gray header outside the panel grid.
- Use sentence case. Avoid unexplained abbreviations inside panels.

## Lines, markers and axes

- Primary data line: 1.2 pt; auxiliary/reference line: 0.7–0.9 pt.
- Axis spine: 0.7 pt; show left and bottom spines only unless a heatmap frame is
  necessary.
- Marker diameter: 4–5 pt for summaries and 2.5–3.5 pt for individual values.
- Error bars: 0.8 pt with 2.5 pt caps.
- No default background grid. A light horizontal guide (`#D9D9D9`, 0.5 pt) is
  allowed only when it materially aids value comparison.
- Start proportion axes at zero unless the panel explicitly displays a narrow
  deviation around a reference and the truncation is conspicuously marked.
- Use log axes for feature counts, sequence-per-RFU distributions and runtime
  only when multiplicative differences are the question.
- Show at most two meaningful decimal places for proportions; use integer
  counts with thousands separators.

## Color system

### Core semantic colors

| Meaning | Color | Hex |
|---|---|---|
| RFU | deep blue | `#0072B2` |
| Threshold-qualified RFU | lighter blue | `#56B4E9` |
| Exact CDR3 | muted orange | `#D89000` |
| Scirpy clonotype | green | `#009E73` |
| TRBV+TRBJ | muted purple | `#8E6C8A` |
| CDR3 length | medium gray | `#7A7A7A` |
| Diversity | brown-gray | `#A67C52` |
| Null/reference | light gray | `#B8B8B8` |

All comparator marks use the same size, opacity and line weight as RFU. RFU is
identified by a stable hue, not emphasis. For monochrome reproduction, combine
color with marker shape and/or line style.

### Dataset colors

| Dataset | Hex |
|---|---|
| Wells | `#4C78A8` |
| GSE190905 | `#F28E2B` |
| GSE157007 | `#59A14F` |
| Scirpy wu2020_3k | `#B279A2` |

Dataset colors are used only when dataset identity is the comparison. Do not
mix dataset and representation color semantics in one legend.

### Continuous scales

- Sequential: `viridis` or `cividis`.
- Diverging: `RdBu_r` centered on an explicitly meaningful reference.
- Never use rainbow, jet or red–green-only scales.
- Heatmap missing values are light neutral gray (`#EFEFEF`), not zero-colored.

## Comparator order

The order is fixed throughout the paper:

1. RFU;
2. exact CDR3;
3. Scirpy clonotype;
4. TRBV+TRBJ;
5. CDR3 length;
6. diversity.

If a representation is unavailable, leave it absent rather than changing the
relative order of the remaining methods.

## Assignment-policy conventions

- Nearest RFU: solid line or filled marker.
- Threshold-qualified RFU: dashed line or open marker.
- Frozen threshold `0.6`: thin dark-gray reference line when relevant.
- Threshold failure is labeled “below frozen threshold,” never “OOD” or
  “unknown.”

## Panel lettering and layout

- Letters sit 2–3 mm left of and 1–2 mm above the plotting area.
- Use a consistent 3–4 mm gutter between adjacent panels.
- Align panel baselines and axis-label edges within each row.
- Schematic panels use the same letter and title placement as plots.
- Avoid more than six panels per main figure; the compressed plan uses four or
  five.
- A panel may contain a small inset only when it answers the same scientific
  question and uses the same source family.

## Legends

- Prefer one shared legend per figure.
- Place legends in unused figure whitespace or beneath the panel row; do not
  cover data.
- Use direct labels when fewer than four series are present.
- Legend order follows visual reading order and the fixed comparator order.
- Dataset, representation and assignment-policy semantics must never share one
  ambiguous color legend.

## Preferred plot forms

| Question | Preferred form |
|---|---|
| Exact parity/invariants | compact matrix or status strip |
| Execution funnel | proportional or aligned boxes with counts; no decorative Sankey |
| Runtime scaling | line/point plot with measured points only |
| Fresh versus cached | paired log-scale bars or dumbbells |
| Feature counts | grouped points/bars on log scale |
| Sparsity versus dimension | labeled scatter |
| Coverage | lollipop or dot plot with zero-based proportion axis |
| Cross-dataset sharing | aligned heatmaps with identical dataset ordering |
| Within versus between | dumbbell or paired summary points |
| Multi-metric comparator | aligned dot plot, identical scale per metric facet |
| Subsampling | line/point curves with fixed fractions |
| Phenotype coupling | abundance-selected heatmap with marginal prevalence |
| Observed versus null | estimate plot with null mean/SD and observed marker |
| Technical lifecycle | pass/fail matrix, not a faux quantitative chart |

Avoid radar charts, pie charts, 3D plots, decorative Sankeys, excessive
gradients and table screenshots.

## Statistical display

- Plot individual biological/sample units when frozen tables contain them and
  the panel remains readable.
- Otherwise display the audited mean/median and dispersion named in the source
  table.
- Do not add significance stars.
- Show empirical permutation probabilities only for prespecified tests and
  label the permutation count.
- Never treat cells, sequences or dependent sample pairs as independent
  biological replicates.

## Export and QA

- Save PDF first, then a 300-dpi PNG preview from the same Matplotlib figure.
- Use `bbox_inches="tight"` only after confirming that panel letters and shared
  legends remain inside the output.
- Set PDF font type to 42 and preserve editable text.
- Verify all text at final physical size, not only when zoomed.
- Check color appearance with a color-vision-deficiency simulator before final
  export.
- Every panel must have a row in `manuscript/figure_source_manifest.tsv` before
  it is considered reviewable.
