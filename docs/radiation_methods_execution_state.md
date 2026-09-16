# Radiation-methods execution checkpoint — 2026-09-16

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
| RP1-14 | Blocked on current asset location | Historical matrix, receptor mappings and donor/visit registry not found. No published-only results fabricated. |
| Wells full standard-backend run | Reused; complete-run and checksum checks passed | 610,429 atlas cells / 24 donors; 303,088 primary-TRB cells / 21 donors; 233,913 threshold cells; 4,928 RFUs. Coverage retains three donors without eligible TRB. |
| GSE190905 cached release | Reused; receptor/RNA/assignment joins verified | 27,655 TCR cells, 43,051 metadata cells, six donors, twelve visits, eight physical library pools. Four source `SBRT`, two `I-SBRT` donors. |
| GSE280982 | Public processed-file and pair inventory completed | Eight tumor visits across three donors have both GEX and TCR entries; three paired blood visits. No external endpoint analysis or new assignment yet. |
| Matos/RfuWAS | Frozen supporting example retained | 17/32 conditional eQTL/caQTL variants; six same-variant overlaps. No formal RFU colocalization. Negative held-out state-prediction result preserved. |

The local asset inventory, full source donor crosswalks and library membership
are external. Historical RP1-14/Wells paths were checked first, followed by
targeted recovery manifests, configs, recovery bundles and document/download
locations. No unrestricted filesystem search or unpublished-cohort analysis was
performed. The user was asked for the current RP1-14 location; it remains pending.

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

Next RP1-14 step requires the current path of the existing RFU matrix and its
published-only sample registry. Inspect actual sample counts, raw/proportion/
transformed values, denominator, assignment policy, reference hashes and
receptor-level availability before specifying a command. A fabricated command
to the absent historical file would not be executable. Do not rerun assignment
to compensate for a missing archive location or reuse mixed-cohort PCA/UMAP.

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

## Validation

Focused multiscale, boundary/library, completed-run and existing longitudinal
tests: **33 passed**. Final full suite: **484 passed, 4 skipped, 32 warnings**
in 22.01 seconds. `ruff check .` passed; `ruff format --check .` reports
**227 files already formatted**. `git diff --check` passed. A repeat invocation
verified the completed cache and returned without recomputation.

The initial full-suite attempt had
483 passes and one missing-manuscript-path failure during the archive/new-draft
transition; that was a transient incomplete documentation state, not a numerical
failure. All 32 dependency warnings and four optional-backend skips are retained.
