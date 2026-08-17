# Manuscript figure and source-table plan

Every panel must be regenerated from a saved source table and manifest. Cells or
sequences are descriptive units, not independent biological replicates.

## Figure 1 — framework and technical benchmark

| Panel | Generating function/workflow | Required source table | Data unit | Uncertainty | Caveat / status |
|---|---|---|---|---|---|
| Workflow/schema | vector assembly from API contract | receptor/result schema tables | field | none | Planned; distinguish original RFU from scRFU. |
| Assignment parity | `examples/original_rfu_parity.py` (planned) | row-level original/scRFU comparison | receptor row | exact mismatch counts | Workflow not yet added; frozen artifacts required. |
| Deterministic chunks | `call_rfu_table` manifests | run/chunk manifest summary | unique query/chunk | repeated-run range if repeated | Serial/parallel synthetic parity implemented; real benchmark pending. |
| Reference coverage | `reference_coverage` | group coverage TSV | sample/cohort/chain | donor/sample summaries | Threshold is not probabilistic confidence or OOD. |
| Runtime/memory | benchmark manifest summary | runtime/RSS/cache table | run | replicate range | 1k/10k/25k bounded runs pending. |
| Downsampling robustness | benchmark utilities | perturbation stability TSV | biological sample/run | seeded replicate distribution | Multinomial abundance resampling is not physical read downsampling. |

## Figure 2 — longitudinal compartment validation

| Panel | Generating function | Source table | Data unit | Uncertainty | Caveat / status |
|---|---|---|---|---|---|
| Design | `validate_longitudinal_design` | design and QC TSV | biological sample | none | Synthetic only; missing visits remain missing. |
| Within/between similarity | `longitudinal_similarity` | tidy pair TSV | sample pair | donor-block interval/bootstrap summary | Do not treat pairs as independent. |
| Donor retrieval | `donor_retrieval` | query-ranking TSV | query sample | donor-level bootstrap | No parameter tuning on evaluation cohort. |
| Persistence/dynamics | `rfu_longitudinal_dynamics` | classifications plus trajectory TSV | donor–RFU–compartment | threshold sensitivity | Descriptive labels, not inferential states. |
| Trajectories | source trajectories | abundance/missingness TSV | biological sample | donor summaries | No hidden imputation. |
| Compartment divergence | `longitudinal_compartment_comparison` | paired comparison TSV | donor-time pair | donor-block bootstrap/restricted permutation | Compartment names are user data, not hard-coded. |

The deep longitudinal cohort is a methods demonstration and cannot support a
general population claim by itself.

## Figure 3 — cross-cohort TCR validation

| Panel | Generating function/workflow | Source table | Data unit | Uncertainty | Caveat / status |
|---|---|---|---|---|---|
| Cohort characteristics | explicit harmonization | mapped fields/QC/missingness | donor/sample | descriptive | Cohorts not yet selected. |
| Frozen coverage | `transfer_cohort` | coverage/score/RFU summaries | sample/cohort | donor/sample summaries | No target refitting. |
| Held-out transfer | held-out manifest plus fixed metrics | registered evaluation TSV | held-out sample | prespecified donor-level method | Blocked until cohort registered. |
| Replication | frozen analysis configuration | cohort-specific effect/source table | donor | model appropriate to cohort | No result claimed in Month 1. |
| Comparator performance | `repertoire_representation`, retrieval/similarity | tidy comparator metric TSV | biological sample/query | donor-block resampling | Comparators are useful baselines, not exhaustive SOTA. |
| Single-cell phenotype interpretation | phenotype coupling/pseudobulk | RFU-by-phenotype/sample TSV | sample and descriptive cell rows | sample-level summaries | Avoid cell-level pseudoreplication. |
| Antigen evidence | VDJdb workflow | pinned-reference evidence/coherence/null TSV | unique sequence/RFU | size-preserving permutations | Evidence is not proof of antigen specificity. |

## Figure 4 — BCR extension (gated and optional)

Only planned after every BCR gate passes. Required source tables would include
heavy/light selection, isotype/SHM/family QC, receptor-specific parity and
coverage, public-cohort results, and held-out validation. Current status:
blocked; not part of the TCR methods claim.

## Optional Figure 5 — unified biological synthesis

Future only. It cannot proceed until independently validated TCR and BCR models
both exist. No generating function or statistical claim is frozen.

## Supplementary source tables

Include API/provenance schemas, all eligibility and exclusion counts, exact
parameters, software/reference hashes, undefined results, comparator settings,
random seeds, environment information, and one panel-to-source-table index.
