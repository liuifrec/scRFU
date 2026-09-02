# Source-data inventory

## 2026-09-01 scverse methods-evidence supplement

The following new source tables remain beside their external public-data
outputs. Paths below are portable labels relative to their evidence roots, not
developer-machine paths. These supplement rows were generated with scRFU 0.4.0.
For the independent public methods project, the historical governed/private
Figure 2 placeholders retained later in this file are archived planning records,
not release or paper requirements; no private cohort was used in this sprint.

| Evidence family | External source table | SHA256 | Analysis unit | Status/caveat |
|---|---|---|---|---|
| Native AIRR scale, 25k | `native_scale/25k/source_table.tsv` | `50873116693ac4660bb5d7ddacedcf496e5cbea69a6f11ca8232ba527c9fc065` | bounded Wells source cell / receptor chain | Complete; X absent; zero table/native and reload mismatches |
| Native AIRR scale, 100k | `native_scale/100k/source_table.tsv` | `ca8ad6a688a43b613e5238302f24d5ce16623f75511bc278460a75566dbb0203` | bounded Wells source cell / receptor chain | Complete; X absent; zero mismatches |
| Native AIRR scale, 250k | `native_scale/250k/source_table.tsv` | `2f65d3de9bc66fb13aeca45a8857a9a6562d16f97a059400dbb36f532f46eded` | bounded Wells source cell / receptor chain | Complete; optional scale gate passed; zero mismatches |
| Genuine Scirpy comparator | `GSE190905/scirpy_comparator/representation_summary.tsv` | `5c052a38c0f393ee76982393b8cefeba40f7f5afa36ee68f05f58df579de7efd` | patient-state sample | Complete; Scirpy 0.22.4 identity/NT, all arms, primary-only dual IR |
| Genuine Scirpy retrieval | `GSE190905/scirpy_comparator/retrieval_summary.tsv` | `1dad4adb237bf4bdc7f8bed8cfd3dc6dbcad9d067c49430a1c92c45f19c6ae65` | held-out patient-state query | Complete; identical candidate sets; descriptive six-donor demonstration |
| Genuine Scirpy downsampling | `GSE190905/scirpy_comparator/downsampling_summary.tsv` | `a09b7e8bdd9c0110a8589f2f97b1e381965faf0e1ae7ee9677353ed3129b2443` | sample representation | Complete; fixed 50%/75% fractions and seeds |
| Representation compression | `representation_consistency/representation_compression.tsv` | `7265a42b43b19992be71f6e9c80c4b8c1e35ac79fb31b799e383568542156b14` | dataset / feature representation | Complete; compression is not claimed as biological superiority |
| Frozen feature sharing | `representation_consistency/cross_dataset_feature_sharing.tsv` | `4306ef66da49c035c3df7109aa03c5d9ad66163dd7ad8a45a5420f3232500878` | dataset pair | Complete; exact CDR3 and RFU spaces compared without outcome fitting |
| Native VDJdb chain linkage | `wu2020_3k/native_vdjdb_linkage/native_chain_vdjdb_summary.tsv.gz` | `f08e54e02a738a37d51a8f07119a0d374bc6734d59fabf7c4eeea01f9f53f54f` | AIRR chain | Complete; annotation evidence only, not antigen specificity |
| Public deeper-longitudinal QC | `GSE345124/qc/sample_qc.tsv` | `acdc40b9c9385498d8a0c6ec0cce6008769f8afb0ae242e88933b3bae7dcf864` | public processed repertoire file | Acquisition/QC complete; RFU execution deferred by bounded gate |

Manifest hashes and public accessions are recorded in
`docs/evidence_index.json`. No H5AD/H5MU, RFU asset, VDJdb database, or public
cohort receptor table is copied into Git.

Inventory freeze: 2026-08-24. Existing runtime files are relative to
`SCRFU_MONTH2_OUTDIR`; this environment variable must point to the external
Month 2 output root. They are not copied into Git. `Pending` means that no
scientific source table exists and the corresponding panel cannot be presented
as evidence.

The software version for all existing rows is scRFU 0.1.0. The RFU reference is
`scrfu-ref-db43f0ffe17d9fe6ec4d0d36231e08b48e52f5b508b7039cc08024fdf0158518`.
Commands are repository entry points with runtime paths supplied through
arguments or environment variables; no machine-specific path is part of the
inventory.

## Main figures and tables

| Item | Source file relative to runtime root | SHA256 / expected checksum | Input dataset | Public/private and shareability | Anonymization | Generation command | Status |
|---|---|---|---|---|---|---|---|
| Fig. 1B original RFU parity | `original_rfu_parity/figure1_original_rfu_parity.tsv` | `9bfc3413abf9d15ae2e1002f2832b2c4e7a97733d8ad776afc5a891808ff2bbf` | adversarial synthetic receptor fixture + official RFU | synthetic/public provenance; shareable | none | C1 below | Complete |
| Fig. 1C input/deduplication | `wells_bounded/wells_sampling_summary.tsv` | `2b2898788c782d0cd15c3bd17481222a8c86b9df1a2d89abca0deac6cab8c197` | Wells public atlas | public-derived aggregate; shareable after source-license review | no participant IDs | C2 below | Complete |
| Fig. 1D–F scaling, parallelism, resume, memory | `wells_bounded/figure1_scaling_source.tsv` | `100bcb6c8dff56256e7e91df3cf7a5598cf30e44717a7c3769e1bcf4e7c276e4` | Wells bounded subsets | public-derived technical table; shareable | none | C2 below | Complete |
| Fig. 1G coverage/policy | `wells_bounded/single_cell_analysis/assignment_policy_summary.tsv` | `f3de40b04831e523b26c1c26700c49a7f9707383464c2125120048321d2ad120` | Wells 25k | public-derived aggregate; shareable after source-license review | aggregate only | C2 below | Complete |
| Fig. 1H robustness | `wells_bounded/figure1_robustness_source.tsv` | `74ea9f21c5c9005a070b9e75b476558aeaf9c342ba09a114d2be8f05e435ceb1` | Wells 25k | public-derived technical table; shareable | none | C2 below | Complete |
| Fig. 2A cohort QC | Pending: `${SCRFU_LONGITUDINAL_OUTDIR}/figure2/cohort_qc.tsv` | Pending; record SHA256 after generation | governed six-volunteer cohort | private-derived; only aggregate public-safe table shareable | suppress source IDs; anonymized analysis IDs only when essential | validation command in `docs/longitudinal_runtime_command.md`, then governed workflow | Blocked |
| Fig. 2B within/between similarity | Pending: `${SCRFU_LONGITUDINAL_OUTDIR}/figure2/within_between_similarity.tsv` | Pending | governed cohort | private-derived; aggregate/pair-class summaries shareable | no donor, sample, visit, or participant identifiers | governed workflow under `docs/methods_freeze.md` | Blocked |
| Fig. 2C donor retrieval | Pending: `${SCRFU_LONGITUDINAL_OUTDIR}/figure2/donor_retrieval.tsv` | Pending | governed cohort | private-derived; aggregate metrics shareable | no query/candidate identifiers | governed workflow under frozen candidate sets | Blocked |
| Fig. 2D persistence/dynamics | Pending: `${SCRFU_LONGITUDINAL_OUTDIR}/figure2/rfu_dynamics_summary.tsv` | Pending | governed cohort | private-derived; category aggregates shareable | no donor-specific RFU trajectory | governed workflow under frozen abundance sensitivity | Blocked |
| Fig. 2E compartment trajectories | Pending: `${SCRFU_LONGITUDINAL_OUTDIR}/figure2/compartment_summary.tsv` | Pending | governed cohort | private-derived; aggregate summaries shareable | generic compartment labels unless disclosure is approved | governed workflow | Blocked |
| Fig. 2F donor/technical robustness | Pending: `${SCRFU_LONGITUDINAL_OUTDIR}/figure2/donor_robustness.tsv` | Pending | governed cohort | private-derived; aggregate resampling output shareable | no bootstrap draw membership or donor ID | governed workflow | Blocked |
| Fig. 3A–C cross-cohort transfer | `cross_cohort/cross_cohort_technical_summary.tsv` | `51d7c0ebb88e4944faabd85b62db7380dd01b7f523d5957d8f5f48fc6a11da2f` | Wells, GSE190905, GSE157007 | public-derived aggregate; shareable after source-license review | none | C3 below | Complete |
| Fig. 3D GSE190905 paired structure | `cross_cohort/gse190905_within_between_comparators.tsv` | `8911ba5abd980c321d17b13dbb5f1cc6028748e5d8c08718131c74d13021b8e5` | GSE190905 | public-derived aggregate; shareable | no source patient label in public-ready view | C3 below | Complete |
| Fig. 3E GSE190905 retrieval/comparators | `cross_cohort/gse190905_donor_retrieval_comparators.tsv` | `d2f79d96b4e0ea39bef170d1fe2a37ff9b5aec95e314f70c993f2d23f6ab468f` | GSE190905 | public-derived aggregate; shareable | no source patient label in public-ready view | C3 below | Complete |
| Fig. 3F held-out robustness | `public_data/GSE157007/heldout_validation/subsampling_stability_summary.tsv` | `b5cbcbb1edc3ee6fb1806092833a66af695bb009391f21f8d87d8d5964f7c866` | GSE157007 | public-derived aggregate; shareable | none | C4 below | Complete |
| Fig. 3G–H Wells phenotype coupling | `wells_bounded/phenotype_coupling_stability.tsv` | `205a4da7606517ad571f7a0ce0f6a356bef42ca480f6a3a87b772742d0c5d385` | Wells 25k | public-derived aggregate; shareable after source-license review | no cell barcode/donor ID | C2 below | Complete |
| Fig. 4 VDJdb coherence | Pending: `${SCRFU_VDJDB_OUTDIR}/figure4/vdjdb_coherence_summary.tsv` | Pending | release-pinned VDJdb + three receptor datasets | public-derived summary shareable; database itself not redistributed | none | C5 below for every frozen dataset/policy/match/ambiguity combination | Blocked |
| Fig. 4 fallback robustness/transfer | existing Fig. 1 robustness and Fig. 3 transfer sources above | use the immutable hashes above | Wells, GSE190905, GSE157007 | public-derived; shareable subject to sources | none | existing commands above | Available fallback |
| Main cohort-characteristics table | Wells sampling table; GSE190905 and GSE157007 run manifests | see `docs/manuscript_claim_audit.md` | three public cohorts | public-derived aggregate; shareable | none | existing bounded/public workflows | Complete for public cohorts |
| Main software/provenance table | `original_rfu_parity/summary.json` and public run manifests | see claim audit | all executed analyses | technical metadata; shareable | none | existing workflows | Complete |

## Generation command registry

All variables below identify external runtime inputs or outputs. C1–C4 are the
commands that generated the audited source families; executing them again is
unnecessary unless an input or implementation changes.

**C1 — official RFU parity**

```bash
python examples/original_rfu_parity.py \
  --rfu-dir "$RFU_DIR" \
  --outdir "$SCRFU_MONTH2_OUTDIR/original_rfu_parity" \
  --score-tolerance 1e-12 --threshold 0.6 --seed 20260824
```

**C2 — bounded Wells technical, robustness, and single-cell sources**

```bash
python examples/wells_month2_validation.py \
  --input "$WELLS_H5AD" \
  --rfu-dir "$RFU_DIR" \
  --outdir "$SCRFU_MONTH2_OUTDIR/wells_bounded" \
  --sizes 1000 10000 25000 --chunk-sizes 1000 5000 --seed 20260824
```

**C3 — harmonized cross-cohort summaries**

```bash
python examples/cross_cohort_validation.py \
  --config "$SCRFU_CROSS_COHORT_CONFIG" \
  --outdir "$SCRFU_MONTH2_OUTDIR/cross_cohort"
```

**C4 — preregistered held-out evaluation**

```bash
python examples/gse157007_heldout.py \
  --input-dir "$GSE157007_INPUT_DIR" \
  --family-soft "$GSE157007_FAMILY_SOFT" \
  --manifest docs/heldout_gse157007_preregistration.json \
  --frozen-reference "$SCRFU_FROZEN_REFERENCE_MANIFEST" \
  --rfu-dir "$RFU_DIR" \
  --outdir "$SCRFU_MONTH2_OUTDIR/public_data/GSE157007/heldout_validation"
```

**C5 — one release-pinned VDJdb analysis cell**

Run this template for each of the three frozen receptor-result inputs and every
prespecified assignment, match, and ambiguity policy; preserve separate output
directories and combine them only after all cells pass validation.

```bash
python examples/vdjdb_antigen_evidence.py \
  --rfu-sequences "$SCRFU_RFU_SEQUENCES" \
  --rfu-rows "$SCRFU_RFU_ROWS" \
  --vdjdb "$VDJDB_PATH" --vdjdb-release "$VDJDB_RELEASE" \
  --expected-sha256 "$VDJDB_SHA256" \
  --outdir "$SCRFU_VDJDB_OUTDIR/$SCRFU_VDJDB_ANALYSIS_CELL" \
  --match-mode "$SCRFU_VDJDB_MATCH_MODE" --v-gene-mode strip_allele \
  --assignment-policy "$SCRFU_ASSIGNMENT_POLICY" \
  --ambiguity-policy "$SCRFU_AMBIGUITY_POLICY" \
  --n-permutations 1000 --random-state 20260824 \
  --save-permutation-values
```

## Public-safe longitudinal derivative contract

No private source table is made public by default. A release candidate may
include only separately generated tables that:

1. contain no participant, sample, visit, library, cell-barcode, or source-file
   identifier;
2. report cohort-level counts, pair-class summaries, aggregate retrieval
   statistics, dynamics-category totals, compartment-level estimates, and
   bootstrap/permutation summaries;
3. use anonymized local analysis IDs only when a plot cannot be interpreted
   without repeated-measures linkage, with the re-identification key retained
   outside Git;
4. suppress or coarsen cells whose combination of fields could disclose a
   participant; and
5. receive a new SHA256 and disclosure review before inclusion.

The private input files, raw result tables, ID mapping, and Figure 2 working
directory remain external. This inventory does not authorize their movement.

## Checksum and release rule

Existing checksums are immutable evidence locators, not promises that the
external runtime directory will always be retained. Before submission, copy
only approved public-safe derived tables into a release-controlled source-data
deposit, recompute SHA256 values, and update this inventory with deposit IDs.
Do not copy RFU assets, VDJdb, H5AD, or private inputs into the repository or
source-data deposit.
