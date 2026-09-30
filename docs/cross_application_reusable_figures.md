# Reusable cross-application quantitative figures

This second asset pack continues from `c9cd3d41f049c04f27e5a9909b5472180d90fa8d`
on `manuscript/radiation-methods`. It uses completed summary tables only. No
workflow cartoon, conceptual schematic or causal arrow is included. The
completed RP1-14 reusable pack is protected and verified without modification.

## Files and recommended placement

The external pack is:

`/home/liuyuchen/data/scRFU_radiation_methods_20260916/results/cross_application_reusable_figures_v1/`

| Key | Quantitative content | Recommended placement |
| --- | --- | --- |
| A | GSE190905 paired receptor/RFU TV in total T, CD4 and CD8; qualified-cell coverage; unique-receptor, saved depth and dominant-receptor sensitivity | Main paper application; poster block |
| B | GSE280982 ordered visit coverage with missing visits; supported receptor/RFU TV; persistent/union and no-shared-receptor/persistent fractions | Main paper external application; poster block |
| C | Seven frozen, separate RP1-14/GSE190905/GSE280982 TV and cancellation strata, retaining units and donor counts | Main paper synthesis; optional poster replacement |
| D | 623 RFU-QTL variants; 17 conditional-eQTL, 32 conditional-caQTL and six both-layer overlaps; six-variant target/class/count matrix | Separate regulatory paper figure; poster block |
| E | Six held-out donor log losses and positive-is-worse deltas in ordinary and purged designs | Supplementary; optional poster/Q&A backup |
| F | Wells donor/tissue coverage, primary and qualified counts, RFU detection breadth across tissues and source cell-type labels | Paper context or supplement; optional poster panel |

Each A–F figure has manuscript and poster exports at
`figures/{manuscript,poster}/cross_application_{A..F}_{manuscript,poster}.{pdf,svg,png}`.
PDF/SVG files contain vector objects; PNGs are 300 dpi. There are 16 export sets
(48 deliverable files): twelve audience figures, two compact components and
two assemblies. Color/grayscale previews and contact sheets are in `preview/`.

The **poster block** is
`figures/poster/cross_application_poster_block.{pdf,svg,png}`. Its A/B/D layout
is 762 × 609.6 mm with text at least 18 pt at native size. The lower-right
rectangle `[left=0.63, bottom=0.055, width=0.35, height=0.365]` in figure
coordinates is blank for BioRender. These coordinates are also in the export
manifest; the exporter verifies that the corresponding PNG pixels are white.
Reducing the physical block also reduces its font sizes. This is a quantitative
block for insertion into an A0 poster, not a completed poster.

The **manuscript assembly** is
`figures/manuscript/cross_application_manuscript_assembly.{pdf,svg,png}`, at
190.5 × 266.7 mm with text at least 8 pt. It combines compact A/B with C.
`figures/panels/cross_application_{A,B}_compact.{pdf,svg,png}` supplies the
individual components. The fuller A/B audience figures retain sampling and
persistence panels. D/F remain separate supporting figures; E should remain
supplementary. The existing RP1-14 main-figure exports are unchanged.

## Reproduction and traceability

Run from the repository root using the already available plotting environment:

```bash
MPLCONFIGDIR=/tmp/scrfu-reusable-mpl \
/home/liuyuchen/data/scRFU_regulatory_20260909/environment/bin/python \
  -m manuscript.scripts.cross_application_reusable_figures \
  --workspace /home/liuyuchen/data/scRFU_radiation_methods_20260916 \
  --regulatory-workspace /home/liuyuchen/data/scRFU_regulatory_20260909 \
  --out /home/liuyuchen/data/scRFU_radiation_methods_20260916/results/cross_application_reusable_figures_v1
```

Append `--verify-only` to check a completed pack. Normal repeat invocation also
verifies and returns without rewriting. Changed code, configuration, captions,
schema, input hashes or plotting software require a new versioned output
directory. Incomplete/nonempty directories are never silently overwritten.

The runner verifies all 231 outputs listed in eight pinned completion
manifests before reading the 24 selected source tables/reports. The pins and
physical export sizes are in
[`cross_application_reusable_figures_v1.json`](../manuscript/config/cross_application_reusable_figures_v1.json).
They include every output in the protected RP1-14 pack. Importing shared
typography, canvas and hashing helpers from its scripts does not run its
analyses or rewrite its figures.

The pack contains 19 plotting TSVs with six per-figure READMEs. Every column is
defined using the checked-in
[`table schema`](../manuscript/figures/cross_application_table_schema.json).
`manifests/figure_index.json` joins each export to exact inputs, source TSVs,
captions and reproduction commands. `manifests/provenance.json` records code,
input and upstream hashes, the software environment and the generating parent
commit. `completion.json` hashes every other output, including the validation
report. Per-figure manifests and READMEs allow traceability without consulting
the manuscript.

Reusable wording is stored in
[`cross_application_reusable_captions.json`](../manuscript/figures/cross_application_reusable_captions.json)
and exported as `manifests/captions.{json,md}`: short titles, one- or two-sentence
poster captions, manuscript legends, take-home messages and placement advice.
Only code, safe aggregate metadata, captions, tests and documentation belong in
Git; no biological tables or rendered bulk outputs are committed.

## Frozen interpretation rules

- GSE190905 uses the cached six-donor release. All displayed total/CD4/CD8
  pairs satisfy the unchanged 100-qualified-cell-per-visit rule. Sensitivity
  draws are summarized within donor; 50 technical draws do not become 50
  donors. D06's frozen common depth is 404 cells, the other five use 500.
- GSE280982 retains all 18 registered donor/tissue/visit combinations, including
  seven missing visits, and all 18 interval combinations with eight supported
  measurements. Missingness is never zero. Visit positions are categorical;
  intervals reuse donors. Persistence uses one-cell detection; the five-cell
  results, including undefined zero-denominator fractions, remain in the TSV.
- Cross-dataset rows preserve sequencing-read, primary-TRB-cell and
  gene-expression-matched primary-TRB-cell units. Cancellation medians are
  copied donor-difference medians, not differences of plotted TV medians.
  No effect is pooled and biological replication is not implied.
- Regulatory matrix rows are six unique variants. Six expression and ten
  chromatin target/rank records are deduplicated across RFU joins. RFU counts
  overlap across variants. TCR and non-TCR expression targets retain their
  source classifications. Overlap is not formal colocalization or causality.
- Prediction uses saved clone-weighted metrics and equal-donor means. Positive
  delta is worse. Ordinary and purged designs are distinct; the purged design
  removes training amino-acid CDR3/V/J identities shared with the test donor.
  Cell-weighted sensitivities are retained. No model or prediction is rerun.
- Wells keeps all 24 atlas donors, including three without primary TRB. Missing
  atlas combinations and observed zero-primary-TRB samples remain distinct.
  Tissue/cell-type breadth is a support summary, not a specificity test. Its
  frozen observed-receptor support field uses amino-acid CDR3 plus V, which is
  explicitly distinguished from the longitudinal nucleotide CDR3/V/J endpoint.

No assignment, RP1-14 endpoint, QTL search or held-out model was rerun. RFU
persistence does not establish preserved antigen function. The paused
144-person cohort remains outside scope.

## Validation

The export run checks source completion and file hashes, frozen counts/medians,
assignment denominators, cancellation identities, pair support, missingness,
persistence fractions, regulatory record deduplication, prediction sign and
saved metric agreement, and Wells support/count consistency. It round-trips all
19 source TSVs, checks every figure reference, rejects raster SVG content,
checks 300-dpi PNG dimensions, audits label bounds/overlap and minimum font
sizes, and checks the reserved poster region.

Focused tests:

```bash
/home/liuyuchen/data/scRFU_regulatory_20260909/environment/bin/python -m pytest \
  tests/test_cross_application_reusable_figures.py \
  tests/test_rp1_14_reusable_figures.py \
  tests/test_radiation_methods.py \
  tests/test_regulatory_rfu_state_validation.py
```

The [safe checkpoint](../manuscript/figures/cross_application_reusable_checkpoint.json)
records final completion/export hashes, independent export reproduction,
vector-PDF checks and no-write repeat verification.
