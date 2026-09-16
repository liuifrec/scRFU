# Radiation-methods execution checkpoint — 2026-09-16

## Current finishing checkpoint

The external application and bounded repair supplement are **complete**. Start
from this section; pending GSE280982 descriptions in the historical checkpoints
below record the earlier state. No RP1-14 endpoint, GSE190905/Wells analysis,
prediction model or QTL discovery was rerun. The generic package is unchanged.

- GSE280982: 13,077 GEX-matched primary-TRB cells; 9,827 qualify at 0.6,
  representing 2,879 RFUs. All 11 released visits meet the unchanged 100-cell
  rule. Seven tumor pairs in three donors and one blood pair are measured;
  ten registered interval combinations lack a required visit. Missing remains
  unavailable, not zero. Radiation-day qualified depths are 113, 143 and 157.
- Tumor RFU TV: pre → radiation-day 0.8897/0.9175/0.8875;
  radiation-day → six weeks 0.9307/0.9194; pre → six weeks 0.7159/0.7308.
  Blood pre → six weeks: 0.6803. All 96 repeated measurement rows contract.
  The exact receptor TV/cancellation values and denominators are in the
  manuscript and `external_donor_results.tsv`; pairs are not independent people.
- At one-cell detection, 67.1%/54.9% of persistent RFUs in the longer tumor
  comparisons lack a shared observed receptor. At five-cell detection, four
  of five radiation-day intervals have no persistent RFUs; the fraction with
  no shared receptor is then undefined. Depth/weighting materially changes
  distances, and matched controls give similar cancellation. No functional
  stability, antigen persistence, absolute depletion or causal effect claimed.
- RP1-14 repair supplement: all 359,530 RFU/score-vector entries unchanged;
  330,433 changed and 29,097 unchanged AA labels. RFU/score-conflicting AA
  groups fall from 41,572/41,581 to zero. Stratified parity agrees for
  **1,600/1,600** RFU labels and threshold decisions, including 485 scores in
  [0.59, 0.61); maximum score error 1.055e-15. This supports the repair, not
  proof of an unavailable historical invocation manifest.
- Bounded author metadata check found no directly reusable external
  barcode-to-state table. The 25-sheet source workbook and published R code
  supply summaries/selected receptors/cluster annotation logic, not the missing
  complete barcode mapping. No RNA atlas was rebuilt.
- Main figures are now framework, frozen RP1-14, existing GSE190905 and external
  GSE280982. Depth controls and repair QC are supplementary. Wells/QTL/negative
  prediction remain supporting evidence, without a forced composite Figure 5.

Active external directories under `$SCRFU_METHODS_DIR`:

| Stage | Status / output hashes | Completion SHA256 |
| --- | --- | --- |
| `prepared/gse280982_rfu_v1` | Complete; 7 top-level outputs plus validated 15 backend chunks | `1c0fbc4d0fdc525719cddfcec0fa209114668de28a0e46a82e47e724488b8a1f` |
| `results/gse280982_v1` | Complete; 15 output hashes | `f285da99bc0ee955d2ca357d75c4c46fd59fcff27f485bfebf7c67cb1c75a75d` |
| `results/rp1_14_repair_qc_v1` | Complete; 16 output hashes, including parity assay tables | `1841c4e6b72d8a870da5e65c0bd034ef2ded17b1a12808dab4b8163b2a72403b` |
| `results/manuscript_finish_v1_1` | Complete; 15 output hashes; final figures and cross-application table | `b0b0e8597b0c7b877c8bcfc25a6f3901ba69d1502dfb5d55a9ddb4243ca1130c` |

Every parent completion/output hash verifies. All four repeat invocations return
without recomputation. Final figures were visually inspected. The first
`manuscript_finish_v1` export is preserved; `v1_1` only repairs diagram readability
and assembles from the same measured source rows. RP1-14, development and
regulatory checkpoint hashes are unchanged. Inputs, config, function/core-code
hashes, executable paths, imported `scrfu` path and package/R versions are in each
provenance file. Analysis ran from parent `d51586f` plus the explicitly hashed
working-tree scripts committed with this finishing checkpoint.

Only two additional public resources were downloaded for the metadata check:
`41467_2025_60827_MOESM4_ESM.xlsx` (68,546,497 bytes; locally computed MD5
`ce08a1c4b1d2978a3f4f80f4f51c147b`) and `author_analysis.R` (85,279 bytes;
published MD5 **verified** `d3f76420491638247349bf0ef13e956d`). Their SHA256s,
source URLs, workbook inventory and availability conclusion are in
`manuscript_finish_v1_1/metadata_availability.json` and provenance. The 22
previously downloaded processed inputs were reused. No QTL, reference or
expression-matrix download was repeated.

### Exact verification / continuation commands

```bash
SCRFU_METHODS_DIR="$HOME/data/scRFU_radiation_methods_20260916"
SCRFU_PY="$HOME/data/scRFU_regulatory_20260909/environment/bin/python"
export MPLCONFIGDIR=/tmp/scrfu-radiation-mpl
export NUMBA_CACHE_DIR=/tmp/scrfu-radiation-numba-cache
"$SCRFU_PY" -m manuscript.scripts.gse280982_analysis assign --workspace "$SCRFU_METHODS_DIR"
"$SCRFU_PY" -m manuscript.scripts.gse280982_analysis measure --workspace "$SCRFU_METHODS_DIR"
"$SCRFU_PY" -m manuscript.scripts.rp1_14_repair_qc --workspace "$SCRFU_METHODS_DIR"
"$SCRFU_PY" -m manuscript.scripts.radiation_methods_finish --workspace "$SCRFU_METHODS_DIR"
```

These commands verify complete runs. Changed scientific inputs/settings/code
fail closed rather than silently reuse or overwrite. Interrupted incomplete
assignments can resume validated backend chunks. No completed stage needs a
new biological run. Current output/source-table indices and manuscript claims
agree; `manuscript/radiation_methods_claims.tsv` distinguishes verified,
descriptive, bounded and unsupported claims.

Validation: **25 focused tests passed** (external preparation/analysis, repair
QC and multiscale); full suite **517 passed, 4 skipped, 32 warnings in 22.15 s**.
`ruff check .` passes; `ruff format --check .`: **238 files already formatted**;
`git diff --check` passes. Log: `logs/finishing_full_pytest.log`. The unchanged
optional skips/dependency warnings are recorded, not suppressed.

Remaining submission work: author/scientific review; authorship, funding and
ethics/access wording; RP1-14 permission/access confirmation and approved
source-table distribution; journal reference/format requirements; a persistent
authorized artifact archive. No new module, positive prediction result, QTL
mediation claim or missing external state analysis is needed to finish this
bounded manuscript. The next informative biological study would add deeper
repeated sampling with actual cell/template counts. Do not retune these results.

## Historical development and recovery checkpoints

Branch: `manuscript/radiation-methods`, created from completed, pushed
regulatory checkpoint `d62f1c1da3891f0553e08e405769c2c8b38b3316`.
The historical `manuscript/scverse-method-figures` branch and its figures are
unchanged. Its `_figure_common.py` helper alone was reused. The former manuscript
and plan are preserved in `manuscript/archive/`; the active specification is
[`radiation_methods_specification.md`](../manuscript/radiation_methods_specification.md).

## Completed and usable

External workspace: `$HOME/data/scRFU_radiation_methods_20260916`.
Active output: `results/development_v1_1/`. The earlier `development_v1/`
remains available but misses three zero-TRB Wells donors in its coverage table;
it is superseded for Wells coverage. All GSE190905 source TSVs are byte-identical
between these versions. No receptor assignments or prediction models were refit.

| Asset | Reuse status | Validated scope |
| --- | --- | --- |
| RP1-14 | Restored, audited and analyzed | Six donors × three visits × CD4/CD8; historical vector alignment repaired without reassignment; 216/216 bounded parity. Six-panel figure and all source tables complete. |
| Wells full standard-backend run | Reused; complete-run and checksum checks passed | 610,429 atlas cells / 24 donors; 303,088 primary-TRB cells / 21 donors; 233,913 threshold cells; 4,928 RFUs. Coverage retains three donors without eligible TRB. |
| GSE190905 cached release | Reused; receptor/RNA/assignment joins verified | 27,655 TCR cells, 43,051 metadata cells, six donors, twelve visits, eight physical library pools. Four source `SBRT`, two `I-SBRT` donors. |
| GSE280982 | Processed inputs acquired and prepared after RP checkpoint | 22 verified contig/barcode files, 13,077 GEX-matched primary-TRB cells; eight tumor and three blood visits. Frozen configuration saved; RFU assignment/endpoints pending. |
| Matos/RfuWAS | Frozen supporting example retained | 17/32 conditional eQTL/caQTL variants; six same-variant overlaps. No formal RFU colocalization. Negative held-out state-prediction result preserved. |

The local asset inventory, full source donor crosswalks and library membership
are external. Historical RP1-14/Wells paths were checked first, followed by
targeted recovery manifests, configs, recovery bundles and document/download
locations. No unrestricted filesystem search or unpublished-cohort analysis was
performed. The initial checkpoint lacked RP1-14 files; the user subsequently
restored both source and RFU directories. That blocker is now resolved.

The reusable public runner and configuration generated two real longitudinal
figures, five indexed panels and source tables for TV decomposition, coverage,
paired similarities, CD4/CD8/state strata, observed RFU persistence, clone
dominance, unique-receptor weighting, depth sensitivity and fixed grouping
controls. Wells tables use all eligible RFUs, with a declared support rule
(at least 50 cells across four donors; 1,471 RFUs meet it).

## Scientific checkpoint

Across six donors, median primary-TRB TV / RFU TV / cancellation is
**0.8399 / 0.7332 / 0.0986** on threshold-assigned cells. TRBV and V/J median
group distances are 0.2419 and 0.4524; differing granularity prevents a claim of
superiority. The inequality holds for all 462 repeated analysis rows.

Fixed feature/size-matched random groupings yield similar cancellation for
several donors. Empirical depth and clone weighting substantially affect
distances. These support measuring aggregation and sampling explicitly, not
equating RFU persistence with functional recovery. The donor, not the cell,
sample pair or resampling replicate, is the biological unit.

The two completed held-out prediction designs are unchanged. Equal-donor,
clone-weighted log loss is 1.850996 → 1.855212 for ordinary donor generalization
and 1.852679 → 1.856996 after excluding shared train/test receptor identities.
Both have six donors and positive deltas (deterioration); neither demonstrates
an RFU prediction gain. See the frozen regulatory report for all fold results.

## Integrity and completion

- Active input/configuration/code/software fingerprint:
  `5b375794f06db1c5dc977295ae36987e43e5268727115f323d4dff3af2b1e717`.
- `completion.json` SHA256:
  `4b5c1fd1f6d9c78e11612547809ce4c806c677835744bf4d51e4cd776171bf17`.
- `evidence_counts.json` SHA256:
  `858906683aeeefed6c347b67805483713ade5809d4aea55e9026af6b05c0da8a`.
- `gse190905_multiscale_pairs.tsv` SHA256:
  `0bf7e912d297457fd6271f4d7239a9a5f6a05e307587e4b178e6c4cf9691866e`.
- `figure_source_index.tsv` SHA256:
  `7ea94f4476c3a9ba9cb807bc4db42ba1c36a98330910eecfd42c052b5ecf5c96`.
- All **32 output hashes** verified. Figures inspected visually; their source
  tables, denominators and limits agree with the manuscript.
- Wells source H5AD's citation identifies dataset version
  `965008a7-e698-413c-a746-8855a945a7c5`, DOI
  `10.1038/s41590-025-02241-4`. The inherited full-file SHA256 is
  `de47634244ff060af7273e4485816326734b36dbabe78253b8c2e70bbf0181bc`.
  Source expression matrices were not reread or reprocessed.

Public source files remain in their existing workspaces. Only small GEO metadata
SOFT files were downloaded. The current runner pins the GSE190905 source tables
and both assignment tables, checks completed-run provenance, and retains the
reference/code/trimer hashes. The persistent regulatory Python environment
imports this repository's `src/scrfu`; executable and package versions are in
`provenance.json`. No `/tmp` dependency environment is required. Generation used
the parent git checkpoint plus explicitly hashed working-tree analysis code;
the validated code is committed with this manuscript checkpoint.

## Commands and unfinished work

Verify/reuse the completed run without rerunning biological measurements:

```bash
SCRFU_DATA="$HOME/data"
SCRFU_METHODS_DIR="$SCRFU_DATA/scRFU_radiation_methods_20260916"
SCRFU_PY="$SCRFU_DATA/scRFU_regulatory_20260909/environment/bin/python"
NUMBA_CACHE_DIR=/tmp/scrfu-radiation-numba-cache "$SCRFU_PY" \
  -m manuscript.scripts.radiation_methods \
  --gse "$SCRFU_DATA/scrfu_public_validation_20260825/GSE190905" \
  --wells "$SCRFU_DATA/wells/scrfu_full_validation_20260825" \
  --sources "$SCRFU_METHODS_DIR/sources" \
  --out "$SCRFU_METHODS_DIR/results/development_v1_1"
```

Expected message: `Completed manuscript stage verified; reused without recomputation.`
Changed scientific inputs/code/settings require a **new** output directory.

RP1-14 preparation and result commands are in the
[reuse audit](rp1_14_reuse_audit.md). Both directories are now complete and their
hashes verify. Repeating those commands returns after verification, without
repeating assignment or endpoints. An explicit `--refresh-exports` on the RP
analysis runner reuses the independently hashed measurement stage; changed
measurement inputs or functions are rejected. Do not use mixed-cohort PCA/UMAP.

The external processed-file inventory lists exact HTTPS URLs and paired visits:

```bash
"$SCRFU_PY" -c 'import pandas as pd; import sys; x = pd.read_csv(sys.argv[1], sep="\t"); print(x[x.assay.eq("TCR") & x.longitudinal_hypr_sample][["accession", "donor", "visit", "compartment", "url"]].to_string(index=False))' \
  "$SCRFU_METHODS_DIR/results/development_v1_1/gse280982_processed_file_inventory.tsv"
```

Acquire only those processed contigs plus matched GEX barcodes needed for the
external application. Reuse the frozen development configuration. RNA state
annotation requires a separate source-label availability check; do not rebuild
a large expression atlas or assume missing patient visits. Existing paired-
receptor inputs genuinely do not yet have RFU assignments for this study.

## Earlier development validation

Focused multiscale, boundary/library, completed-run and existing longitudinal
tests: **33 passed**. Final full suite: **484 passed, 4 skipped, 32 warnings**
in 22.01 seconds. `ruff check .` passed; `ruff format --check .` reports
**227 files already formatted**. `git diff --check` passed. A repeat invocation
verified the completed cache and returned without recomputation.

The initial full-suite attempt had
483 passes and one missing-manuscript-path failure during the archive/new-draft
transition; that was a transient incomplete documentation state, not a numerical
failure. All 32 dependency warnings and four optional-backend skips are retained.

## Completed RP1-14 checkpoint

The restored raw directories remain in the checkout, untracked and protected by
exact ignore rules. `git ls-files data` lists only `data/schema.md`. No biological
input is staged. The detailed external inventory has **81 file/member records**
(45 physical files plus 36 ZIP members). Source crosswalks remain external.

Active RP1-14 outputs:

- `prepared/rp1_14_v1/`: raw-count reconciliation, corrected source-row mapping,
  sample registry, reference parity and preparation completion manifest.
- `sources/rp1_14_asset_manifest.tsv` and `sources/rp1_14_provenance.json`: detailed
  inventory and supplemental audit of the recovered map-aware backend. R 4.1.2
  reports its version to stderr; the supplemental provenance captures it explicitly.
- `results/rp1_14_v1/`: measured distances, coverage, similarities, RFU persistence,
  dominant-clone/read-depth/fixed-map sensitivities, six-panel PDF/PNG,
  `figure_source_index.tsv`, provenance and final `completion.json`.

**Verified findings:** six donors, 36 samples, 14.95–24.86 years between first/last
visits; 4,342,203 raw sequence rows and 70,960,689 reads. The historical matrix
is assignment-row mass per 10,000, not read abundance. Its map-aware export
misaligned 330,433 sequence labels; the saved RFU/score vectors were preserved
and realigned with source rows. All repeated-sequence conflicts disappeared;
216/216 parity receptors agree with the current reference.

The fixed reused dictionary covers 617,506 qualified sample–receptor observations,
508,972 distinct nucleotide-CDR3/V/J identities and 4,958 RFUs. Qualified reads are
28.35–78.44% of productive source reads. Earliest/latest median receptor/RFU TV is
0.8246/0.3055 in CD4 and 0.7180/0.4309 in CD8. All **432** repeated measurement
rows satisfy contraction. CD8 RFU TV is higher in five of six paired donors;
median paired difference +0.1192. These are descriptive six-person results.

Matched-map cancellation is similar to the observed grouping. At five-read RFU
detection, 41.4%/59.8% of persistent RFUs have no shared observed receptor
(median across CD4/CD8 donors). Aggressive 500-read subsampling and unique-clone
weighting materially change distances. No RFU-specific functional stability,
antigen persistence or population aging model is inferred.

The original `AssignRFUs_with_map` function was recovered in the local historical
RFU checkout and contains the exact truncation error. Its hash and bounded parity
support reuse; the original invocation/reference execution manifest remains
unavailable. No full RFU assignment was rerun. No 144-person inputs, aging-score
shortlists or mixed embeddings enter the RP pipeline.

Integrity:

- Completion SHA256: `6917164a67eca9420201ff78ac9ab7e889fe342bc25fddc896a9130628aa8a56`.
- Evidence counts: `8c33c0e9d56feb59d8b04965ee8f029bce262c5a1e9364326038402031f15331`.
- Main distance table: `e30ef684fbab8c1b3f75a1458ab84dd29e6c8dc7543447fc87ed9350af109dd7`.
- Figure source index: `342e9c7fd588a45f747e5ba1d0fb3d27ada95ff6cf8054a0a56ea71c6bb20e55`.
- Figure PDF: `45d87bacf433583591f7179078c5aa49593b0abbcc107a6195adc26cbb362808`.

Two export-only errors during execution (canonical score-column name, then a
filtered aggregate count field) were corrected. Completed preparation and
measurement stages were verified and resumed; biological distances, controls
and the bounded reference assay were not repeated. Final figures were inspected
visually after moving legends away from data. The manuscript, claims and source
index agree. The existing GSE190905, Wells, prediction and QTL results are unchanged.

Current validation: **37 focused tests passed**; full suite **503 passed,
4 skipped, 32 warnings** in 22.63 seconds. `ruff check .` passes;
`ruff format --check .` reports **231 files already formatted**;
`git diff --check` passes. Logs are external under `logs/rp1_14_*`.

Next: checkpoint this completed RP figure, then acquire the already inventoried
GSE280982 processed contigs and matched GEX barcode metadata. Verify actual
productive receptor/cell and donor/visit coverage before assignment. Keep the
frozen GSE190905 configuration, missing visits and blood/tumor denominators.

## GSE280982 source checkpoint after RP1-14

RP1-14 checkpoint **`7aa878285ccba53f44367521a24a20b7b4ac9830`** was committed and
pushed before starting the external source preparation. Original public
GSE190905/Wells/QTL/prediction outputs remain unchanged.

All **22** requested processed files are now available (11 filtered contig tables
and 11 matched GEX barcode lists), totaling **3,094,261 bytes**. Gzip CRC checks
pass and local SHA256 values are recorded. The inventory did not provide a
publisher checksum; none is invented. No expression matrices, raw sequencing
reads or unrelated cohorts were downloaded.

`prepared/gse280982_v1/` contains a verified sample registry, canonical primary-TRB
receptors, source manifest, original development configuration snapshot and
completion manifest. It represents **13,077 cells / 7,086 unique CDR3 amino-acid
sequences**, with nucleotide CDR3 and V/J present for every selected receptor.
All selected primary-TRB barcodes match the corresponding GEX sample. Multiple-
TRB cells are recorded; the existing adapter chooses primary chains by productive/
high-confidence status, UMI count, reads and source order. The original adapter
code hash is pinned. Barcode IDs are namespaced by sample before concatenation.

There are **11,418 tumor cells across eight visits and three donors**, and
**1,659 blood cells across three visits and two donors**. Only one blood donor
has two available visits. Tumor visit counts are 3/2/3; one donor's final
visit is missing. The first published donor has GEX entries without matching
released TCR entries in the pinned inventory; its TCR counts remain UNKNOWN.
The three observed last-radiation-day tumor samples have **170, 185 and 215**
primary-TRB cells before RFU qualification. Do not lower the frozen 100-qualified-
cell rule if post-assignment counts fall below it.

Source sample treatment strings independently establish pre-treatment, last day
of radiation and six weeks post-radiation; one `Radition` typo is normalized
explicitly. The [publication](https://doi.org/10.1038/s41467-025-60827-w) confirms
the treatment/biopsy order. Cell counts here measure processed receptor coverage,
not absolute lymphocyte depletion. No external RFU endpoint has been examined.
The downloaded inputs do not include author cell-state labels; obtaining a
compatible annotated metadata table remains separate from receptor assignment.

Hashes:

- Preparation completion: `8f346bc4e4363b8d8d99960d54119854c8851c57be84111b514990427e3112cc`.
- Source manifest: `9a569c6963f197f3242a82555bbdc1127525638b5b69a5570b9ef44a50d8eaa3`.
- Verified registry: `8443d723405eace203ab1fb95349e6a0089bae9d9a0ab8f9892c14088e1a24e1`.
- Canonical receptors: `5303df577df6256297870126093d6929a33b6d086a088b4e3bada7fe03d3dbd5`.
- Frozen configuration snapshot: `cb51b5933bb1494401934f016f8f964ef4327e62f842d3d931d575587fa89e8b`.

Revalidate without downloading or rebuilding:

```bash
NUMBA_CACHE_DIR=/tmp/scrfu-radiation-numba-cache MPLCONFIGDIR=/tmp/scrfu-radiation-mpl "$SCRFU_PY" \
  -m manuscript.scripts.gse280982_prepare \
  --inventory "$SCRFU_METHODS_DIR/results/development_v1_1/gse280982_processed_file_inventory.tsv" \
  --pairs "$SCRFU_METHODS_DIR/results/development_v1_1/gse280982_public_pair_availability.tsv" \
  --raw "$SCRFU_METHODS_DIR/sources/GSE280982" \
  --out "$SCRFU_METHODS_DIR/prepared/gse280982_v1"
```

Use `--download` only when the inventoried inputs are absent. Downloads use two
workers, bounded retries and resumable partial files. This checkpoint's repeat
invocation verified completion and returned without repeating preparation.

**Exact next analysis:** load `gse280982_primary_trb_matched_gex.parquet` and call
`scrfu.tl.assign_rfu` with the pinned standard reference, threshold 0.6,
deduplication, chunk size 500 and an external resumable work directory. Preserve
primary-chain choices and sample labels. Validate completion/coverage, then apply
the frozen 100-cell minimum and multiscale measurements separately to tumor and
blood, with missing visits unchanged. Independent RNA labels are needed for
CD4/CD8/state-stratified endpoints; no new classifier search is planned.

Latest validation (including this preparation code): **4 new focused tests
passed; full suite 507 passed, 4 skipped, 32 warnings** in 21.95 seconds.
`ruff check .` passes; `ruff format --check .` reports **233 files already
formatted**. RP-only validation above remains its original 503-test checkpoint.
