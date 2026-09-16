# Radiation-methods application: frozen development specification

Working title: **scRFU: reference-anchored analysis of longitudinal T-cell
repertoire remodeling with applications to radiotherapy**.

This specification and `config/radiation_methods_v1.json` were written before
computing the new longitudinal endpoints. GSE190905 and Wells have already
been used in package development; neither is an untouched validation cohort.
GSE280982 metadata and its publication were inspected before external endpoints.
No external endpoint results were available at the original specification
checkpoint. The finishing application below records execution of these settings;
the configuration itself remains byte-for-byte unchanged.

## Authorized inputs and reuse

| Dataset | Role | Permitted reuse and boundary |
| --- | --- | --- |
| RP1-14 | Published longitudinal CD4/CD8 demonstration | Reuse assignment and sample-wise normalization only after recovering the matrix, actual visit metadata and explicit published/authorized sample membership. No mixed-cohort embeddings or transformed scores. |
| Wells | Public tissue/cellular context | Reuse the complete standard-backend assignment run; all eligible RFUs with a declared public support rule. Do not use the missing historical map-aware survivor-selected subset. |
| GSE190905 | Public development application | Cached six-donor release; reconstruct donor/visit/library linkage from source metadata. Do not impute the seventh donor in the final publication. |
| GSE280982 | Candidate external tumor application | Separate blood/tumor denominators and exact available processed GEX/TCR pairs. Missing visits and unavailable TCR files remain missing. |
| Matos/RfuWAS | Frozen supporting example | Reuse regulatory checkpoint `d62f1c1`; no model refits or genetic searches. Numeric RFU cross-links require reference identity. |

The paused 144-person RERF cohort, its outcomes, RFU shortlists, learned
transformations and mixed-cohort embeddings are excluded. All biological files
and identifying metadata remain outside Git. Published does not automatically
mean redistributable. Pinned input hashes, explicit dataset/sample registries
and authorization checks prevent unregistered inputs from entering this runner.

## Receptor identity and estimands

The development dataset supplies one selected primary TRB per cell. The primary
observed receptor identity is nucleotide CDR3 + source V + source J. It is **not**
a full paired alpha/beta clonotype. The source primary-chain choice is retained;
secondary chains are not silently pooled. RFU assignment remains based on the
cached amino-acid CDR3 result, with original `RFU1`–`RFU5000` labels intact.

Use threshold-qualified cells at the existing 0.6 threshold as the primary
analysis. Compare RFU, TRBV and V/J groupings on those **same cells**. Report
both TCR coverage among metadata cells and threshold coverage among TCR cells.
Nearest assignments on all receptor-bearing cells are a separate sensitivity.
Within each policy, normalize observed counts separately in each visit. Missing
visits are never empty biological samples or imputed zeroes.

For a fixed receptor-to-group map, report receptor total variation, group total
variation and their difference:

\[
D_c=\tfrac12\sum_c|p_c(1)-p_c(0)|,\quad
D_g=\tfrac12\sum_g|\sum_{c\in g}p_c(1)-\sum_{c\in g}p_c(0)|.
\]

The nonnegative difference is **aggregation cancellation**. Contraction is a
mathematical property of aggregation, not evidence of maintained function,
antigen recognition or clinical recovery. Counts or masses are accepted; no
integer counts are inferred from proportions or transformed matrices.

## Development analyses fixed before endpoints

- Biological unit: donor. Report every paired donor before aggregate summaries;
  no independent-sample tests on cells, sample pairs or resampling replicates.
- Pair support: at least 100 retained primary-TRB cells in each visit. Report
  failed strata explicitly. Whole repertoire, broad CD4/CD8 source categories
  and individual source expression clusters are separate strata.
- Primary weighting: cells. Sensitivity: one observation per unique primary-TRB
  receptor per sample. This addresses expansion; it does not estimate unseen
  receptor diversity.
- Depth sensitivity: 50 independently seeded subsamples of actual observed cells
  without replacement; equal pre/post depth per donor, capped at 500 cells.
  These replicate intervals describe subsampling instability, not donor-level
  uncertainty or a sampling-only causal model.
- Expansion sensitivity: remove the union of each visit's most abundant primary
  TRB identity, apply the same exclusion at both visits, renormalize and report
  retained mass.
- Fixed random grouping controls: 30 maps, seeds starting at 20260916. Permute
  RFU group labels among the pooled eligible unique receptor universe within
  exact source V, J and five-amino-acid length bins. This preserves global group
  size and group-by-feature counts, but not per-visit group sizes. Report the
  exchangeable fraction and changed-label fraction. No calibrated null p-value.
- Persistence: observed RFU detection at at least one cell and at least five
  cells; report one-cell fraction `1/N` for each visit and within-RFU observed
  shared receptors, unique receptor counts and dominant-receptor fraction.
  Nondetection does not establish biological absence.
- State composition is restricted to cells with source cluster labels in the
  processed TCR table. Do not impute labels to RNA-only cells.
- Wells support: at least 50 threshold-qualified cells across at least four
  donors, chosen without survivor outcomes. Save all RFU support counts and
  donor/tissue/cell-type tables, including unsupported RFUs for transparency.
- Retain the previously completed negative cell-state prediction benchmark.
  Do not retune it for this manuscript.

## Figure and source contract

| Planned figure | Evidence and status |
| --- | --- |
| 1. Framework | Measurement/evidence schema, with provenance, coverage and interpretation boundaries. No new computational-performance claim. |
| 2. Published RP1-14 trajectories | Frozen six-panel benchmark, with matched controls retained. Repair-QC table is supplementary. |
| 3. GSE190905 radiotherapy trajectories | Existing paired-donor figure; depth/grouping sensitivity figure retained as supplement. |
| 4. GSE280982 external application | Completed donor/tissue/interval measurements; missing visits and low depth explicit. Cross-application source table does not pool effects. |
| Supporting cellular / regulatory evidence | Wells tables, frozen QTL and negative prediction evidence remain supporting material rather than an overcrowded combined Figure 5. |

Generated figures are external. Each panel records source hashes, combined
input/configuration/code fingerprint, command, denominator, biological unit,
interpretation limits and completion status in `figure_source_index.tsv`.
`completion.json` is written last and verifies every output. A different
scientific fingerprint requires a new output directory.

## Restored RP1-14 adaptation, before endpoint inspection

`config/rp1_14_v1.json` records the bulk-read adaptation of the frozen controls.
The source dataset has six donors and three paired-compartment visits. Visits
are sorted by measured collection age; filename suffixes reverse time. Historical
maps require the demonstrated C-start filtering/prefix-alignment repair, never
new RFU values. A bounded 216-receptor reference comparison gates reuse.

The fixed union of existing productive amino-acid assignments is applied to
full productive source counts in every visit, preventing the historical per-visit
top-10,000 selection from creating false absence. The estimand is conditional on
this recoverable universe, with per-sample excluded read mass. Primary pairs are
earliest/latest available visits; all observed pairs are a separate repeated-
measurement sensitivity. Source integer reads replace cell counts; all seeds,
30 grouping controls, 50 subsamples and the 500-observation cap are unchanged.
Grouping uses V/J and five-amino-acid length bins, with unresolved source calls
explicitly marked. No constraints are retuned after seeing RP1-14 endpoints.

These are observed read/receptor summaries, not absolute cell abundance, template
sampling or an aging prediction model. No legacy eRFU score or mixed-cohort
transformation is reused. Historical row-normalized profiles and read-weighted
measurements have different denominators and are never silently interchanged.

## Finishing execution, with no endpoint retuning

The prepared GSE280982 configuration SHA256 is
`cb51b5933bb1494401934f016f8f964ef4327e62f842d3d931d575587fa89e8b`.
External assignment and endpoint stages require this complete snapshot to match
the development configuration. Tumor and blood are separate. All three named
intervals remain separate; required visits must pass the 100-qualified-cell
rule even for nearest-policy sensitivity. No unavailable visit is inferred.
The eight supported pairs and all failed availability combinations are recorded.
Metadata inspection preceded endpoints, so this is not described as untouched
validation. No new RNA atlas or state classifier was constructed.

The RP1-14 supplement checks existing vectors and a fixed 1,600-receptor
stratified parity sample; it changes no biological endpoint. The source
sampling design is documented in `rp1_14_repair_qc.py` and its external
`parity_design.json`. Reporting-only assembly verifies completed parent
manifests before summarizing or drawing panels. The initial and final reporting
exports remain versioned; no model or biological measurement was rerun to fix
plot labels. Exact command and integrity details are in the execution state.
