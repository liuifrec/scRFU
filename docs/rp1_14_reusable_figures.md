# Reusable RP1-14 figures for APBC 2026 and the methods manuscript

This figure-only pipeline reuses the completed RP1-14 results at checkpoint
`9078dc3e9e6e7adb95e7cc675bc5cd5e014f242b`. It does not import the biological
analysis runners, read receptor sequences, rerun assignments or generate new
random maps/read draws. The paused 144-person cohort is outside its scope.

## Deliverables and recommended placement

The external output root is
`$SCRFU_METHODS_DIR/results/rp1_14_reusable_figures_v1/`.
Every A–F figure has a manuscript and poster version in vector PDF, outlined-text
SVG and 300-dpi PNG. The same drawing functions build the assemblies.

| Figure | Content | Poster placement | Manuscript placement |
| --- | --- | --- | --- |
| A | Six-donor/three-visit design, complete CD4/CD8 coverage, donor read depths and global qualified totals | Overview | Main |
| B | Paired receptor/RFU/TRBV/VJ total variation; poster focuses on receptor versus RFU | Alternative main comparison | Main |
| C | Observed aggregation cancellation and 30 fixed matched-map controls for every donor/compartment | Central message | Main |
| D | Persistent RFUs without a shared observed receptor; explicitly hypothetical vector schematic | Central message | Main |
| E | Read versus unique-receptor weights and frozen common 500-read sensitivity | Small inset | Supplement |
| F | Unchanged assignment/score vectors, eliminated conflicts, 1,600-receptor parity and threshold coverage | Optional Q&A backup | Supplement only |

Key artifacts relative to the output root:

- `figures/poster/rp1_14_poster_block.{pdf,svg,png}`: A/C/D/E block,
  **762 × 508 mm**; not a complete poster.
- `figures/manuscript/rp1_14_manuscript_main.{pdf,svg,png}`: refined A–D
  main-figure draft, **190.5 × 264.16 mm**.
- `figures/{manuscript,poster}/rp1_14_{A..F}_{manuscript,poster}.{pdf,svg,png}`:
  twelve standalone figures.
- `figures/panels/`: separate compartment keys and explanatory turnover
  schematics for both audiences. The schematic uses arbitrary a/b/c/d labels;
  it is not a measured example or evidence of antigen recognition.
- `tables/`: fifteen source TSVs and A–F README files defining every column,
  numerator, denominator, input hash and generating command.
- `manifests/`: per-figure metadata, complete source/output indices, provenance,
  export dimensions, font/layout checks, and caption drafts in JSON and Markdown.
- `preview/`: contact sheets plus color and grayscale previews.

The original six-panel `results/rp1_14_v1/rp1_14_longitudinal.{pdf,png}` remains
byte-for-byte unchanged. Its within/between-donor similarity and separate paired
CD4/CD8 RFU-TV views remain available. The new main-figure draft makes
representation, matched controls and observed turnover explicit without changing
the underlying endpoints. The frozen manuscript figure is preserved alongside
this improved presentation option; no existing result directory is overwritten.

## Reproduction

From the repository root, using the existing environment:

```bash
SCRFU_METHODS_DIR="$HOME/data/scRFU_radiation_methods_20260916"
SCRFU_PY="$HOME/data/scRFU_regulatory_20260909/environment/bin/python"
export MPLCONFIGDIR=/tmp/scrfu-reusable-mpl
"$SCRFU_PY" -m manuscript.scripts.rp1_14_reusable_figures \
  --workspace "$SCRFU_METHODS_DIR"
```

Reinvoke with `--verify-only` to check every input/output hash and every saved
source table against its figure-specific projection. A normal repeat also
verifies the complete export without rewriting any files. A changed input,
plotting script, caption file or software version requires a new output path
via `--out`; complete exports cannot be silently replaced. Figure selection
settings are immutable in this v1 runner. Rendering a new layout version does
not require rerunning any biological workflow.

The input lock is
[`rp1_14_reusable_figures_v1.json`](../manuscript/config/rp1_14_reusable_figures_v1.json).
It pins these upstream completion SHA256s and verifies all 35 listed outputs,
including the preserved original figures:

- `rp1_14_v1/completion.json`:
  `6917164a67eca9420201ff78ac9ab7e889fe342bc25fddc896a9130628aa8a56`.
- `rp1_14_repair_qc_v1/completion.json`:
  `1841c4e6b72d8a870da5e65c0bd034ef2ded17b1a12808dab4b8163b2a72403b`.

The runner is
[`rp1_14_reusable_figures.py`](../manuscript/scripts/rp1_14_reusable_figures.py);
all drawing functions are in
[`_rp1_14_figure_drawing.py`](../manuscript/scripts/_rp1_14_figure_drawing.py).
The committed, safe
[caption source](../manuscript/figures/rp1_14_reusable_captions.json)
contains a title, two-sentence-or-shorter poster caption, fuller manuscript
legend, take-home message and placement recommendation for every figure and
assembly. Generated captions are copied to
`manifests/rp1_14_reusable_figure_captions.{json,md}`.

## Scientific and rendering contracts

Primary comparisons remain earliest-to-latest visits, threshold 0.6, and
sequencing-read weights on the frozen qualified reuse universe. Distinct
nucleotide-CDR3/V/J identities are **global unique counts**, not summed sample
richness. Persistence remains at least five qualified reads at both endpoints;
the no-shared-receptor fraction is divided by persistent RFUs, not the detected
union. Undefined ratios remain undefined.

Cancellation is receptor TV minus group TV for each donor, not the difference
between compartment medians. The control ranges are linear 2.5th–97.5th
percentiles of 30 frozen maps, which retain pooled label counts within recorded
V/J and CDR3-length strata. The 500-read sensitivity first summarizes the 50
already-completed draws **within each donor**, then summarizes six donors.
Neither control maps nor read draws are independent people or donor confidence
intervals. The only matched-depth sensitivity shown is the frozen common
500-read-per-visit result; no additional depth analysis was introduced.

Lower group TV partly reflects aggregation. Matched controls produce similar
cancellation; persistence cannot establish preserved antigen function.
Sampling/weighting materially affect the estimates. The repair assay is bounded
computational evidence, not recovery of the missing historical invocation
manifest. RFU means receptor functional unit; TV means total variation; CDR3
means complementarity-determining region 3; V/J are variable/joining calls and
TRBV is the T-cell receptor beta variable-call category.

Fonts are consistent DejaVu Sans, at least **18 pt** in poster exports and
**8 pt** in manuscript exports at native size. The poster block should remain
at least 762 mm wide to preserve its 18-pt minimum; use individual panels for
other layouts. Filled blue circles and open orange squares distinguish CD4/CD8
even in grayscale, with matching bar hatches and line patterns where needed.
PDFs contain vector paths and embedded fonts; SVG text is converted to paths
for portable typography. Each export is checked for text outside the canvas
and overlapping text boxes; color/grayscale previews support visual review.

## Validation

The build verifies all upstream hashes, the 36-sample design, frozen primary
medians (18 checks), shared receptor-TV denominators across representations,
the cancellation identity, both persistence fractions, exact replicate counts,
repair counts and all source-table numeric round trips. It checks every manifest
reference, output hash, PNG pixel dimensions, SVG vector content, font minimum
and text layout. No donor-level data or bulk exports are committed to Git.

Focused tests cover endpoint/policy/weight filters, duplicate/missing pairs and
replicates, within-donor versus pooled technical summaries, immutable selection,
hash tampering (including unused upstream artifacts) and manifest path boundaries.

Completed validation: **28 focused tests passed**; repository-wide Ruff lint and
format checks passed (242 formatted files); `git diff --check` passed. An
independent clean export reproduced all **54 PDF/SVG/PNG files and 15 source TSVs
byte-for-byte**. All 18 PDFs contain one page and no raster images. A repeated
normal invocation left all 130 output-file modification times unchanged.
Final color/grayscale previews and an independent Poppler render of the poster
PDF were visually inspected. The output completion SHA256 is
`0b9cbd974cabd4ea9b84a3852e524fdd0f1338645a8f065391df76f89b3a0d99`;
it protects 129 files. The small, non-biological
[checkpoint record](../manuscript/figures/rp1_14_reusable_checkpoint.json)
also records the six assembly export hashes and the rendering fingerprint.

```bash
"$SCRFU_PY" -m pytest tests/test_rp1_14_reusable_figures.py \
  tests/test_rp1_14.py tests/test_rp1_14_repair_qc.py
"$SCRFU_PY" -m ruff check .
git diff --check
```
