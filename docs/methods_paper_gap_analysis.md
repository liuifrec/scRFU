# Methods-paper gap analysis

Status reflects repository evidence on 2026-08-17. “Synthetic” is not real-data
validation. “Bounded public” refers only to previously validated bounded Wells
workflows, not biological inference. Priorities are P0 (submission-critical), P1
(important), P2 (gated/deferred). Acceptance always includes offline synthetic
tests and documented provenance.

## Core interfaces

| Capability | Status | Public API / implementation / tests / docs | Scientific validation, data requirement, risk, priority, acceptance |
|---|---|---|---|
| Python-native workflow | Implemented | `scrfu.pp/adapters/io/tl/pl`; package tests; API docs | Bounded public receptor workflow; low risk; P0; clean-wheel workflow passes. |
| CLI | Implemented | `scrfu` / `cli.py`; CLI tests; README | Synthetic only for new methods; medium; P0; CLI parity for stable workflows. |
| AIRR input | Implemented | `prepare_receptors`; `adapters.py`; adapter tests/docs | Synthetic; public AIRR cohort needed; low; P0; field/chain parity. |
| AnnData | Implemented | adapter plus `call_rfu`; generic tests; API contract | Bounded public; medium; P0; no expression materialization where selective. |
| MuData | Partially implemented | explicit `modality` routing in `prepare_receptors`; methods tests | Synthetic duck-typed; installed-MuData test needed; medium; P1; modality/alignment fixture passes. |
| scirpy | Implemented | `scirpy_airr`; adapter tests/docs | Synthetic; public object needed; low; P0; supported DataFrame AIRR contract. |
| Bulk receptor tables | Implemented | generic AIRR/DataFrame adapters and table RFU API | Synthetic; public bulk cohort needed; low; P0; sample metadata remains separate. |
| Single-cell receptor tables | Implemented | Cell Ranger, AIRR, Wells; adapter tests | Bounded public; low; P0; stable cell reconstruction. |
| Cell Ranger VDJ | Implemented | `cellranger_vdj`; focused tests/docs | Synthetic; low; P0; CSV/DataFrame/selection parity. |
| Multiple receptor chains | Partially implemented | schema/adapters retain chains | RFU model is reference-specific/TRB; high; P0; reject incompatible model/chain combinations. |
| Sample/donor/tissue/state/compartment/timepoint aggregation | Partially implemented | pseudobulk, phenotype, longitudinal APIs | Synthetic except bounded sample/phenotype paths; medium; P0; explicit keys and no pseudo-samples. |

## Assignment engine

| Capability | Status | Public API / implementation / tests / docs | Scientific validation, data requirement, risk, priority, acceptance |
|---|---|---|---|
| Original-RFU parity | Partially implemented | integration tests and RFU table engine | Exact 1k bounded parity reported; full frozen comparison table pending; high; P0. |
| Exact score/threshold parity | Tested on bounded public data | integration/parity assertions | Extend across edge cases and versions; high; P0; zero unexplained mismatches. |
| Deterministic deduplication/chunking/order invariance | Implemented | `tl_rfu_repo.py`, chunk engine, chunk tests | Synthetic and bounded; medium; P0; byte/row-equivalent outputs. |
| Serial restartability | Implemented | chunk manifests/cache validation | Synthetic and bounded; low; P0; invalid cache never reused. |
| Parallel execution | Tested only synthetically | `max_workers`, process/thread executor; fake-runner tests | Real scaling not run; medium; P1; serial/two-worker/chunk-size equality. |
| Memory efficiency | Tested on bounded public data | selective H5AD/cache readers | Larger bounded scaling pending; high; P0; no `X`/`raw/X` materialization. |
| Multiple-chain behavior | Partially implemented | chain filters and schema | Reference compatibility enforcement incomplete; high; P0. |
| Assignment confidence | Partially implemented | score, threshold, status aliases | Score is similarity, not calibrated probability; high; P0; terminology frozen. |
| Reference coverage / low similarity | Tested only synthetically | `reference_coverage`; diagnostics tests | Real distributions pending; medium; P0; group counts and quantiles reconcile. |
| OOD terminology | Implemented documentation constraint | API freeze and diagnostics docs | No OOD model; high; P0; use “below threshold/low reference similarity.” |
| Frozen-reference provenance | Tested only synthetically | `FrozenRFUReference`; transfer tests | Actual artifact manifest pending; high; P0; hashes and identifier verify. |
| Cell-level reconstruction | Implemented | table engine/backend tests | Bounded public; high; P0; row count/IDs invariant. |

## Analysis

| Capability | Status | Public API / implementation / tests / docs | Scientific validation, data requirement, risk, priority, acceptance |
|---|---|---|---|
| Conventional repertoire metrics | Implemented | `repertoire_metrics`; hand tests/docs | Synthetic; public comparisons pending; low; P0. |
| RFU metrics, pseudobulk, overlap, phenotype coupling | Implemented | `scrfu.tl`; downstream tests/docs | Bounded workflow plus synthetic formulas; medium; P0. |
| Antigen evidence | Implemented | offline VDJdb APIs/tests/docs | Synthetic; pinned local real reference smoke pending; medium; P1. |
| Longitudinal analysis | Tested only synthetically | longitudinal module/tests | No cohort run; high; P0; prespecified deep and public validation. |
| Donor retrieval | Tested only synthetically | `donor_retrieval` | No optimization/evaluation cohort; high; P0; leave-timepoint-out public result. |
| Persistence/dynamics | Tested only synthetically | explicit trajectory classifier | Defaults not biologically frozen; high; P0; sensitivity and source trajectories. |
| Cross-cohort transfer | Tested only synthetically | frozen reference and `transfer_cohort` | No cohort transfer yet; high; P0; development/heldout separation. |
| Held-out validation | Partially implemented | held-out manifest contract | Cohort unregistered; high; P0; manifest predates analysis. |
| Source-table export | Implemented | result tables/cache/workflows | Synthetic/bounded; low; P0; each figure panel traces to table. |

## Receptor modes

| Capability | Status | Public API / implementation / tests / docs | Scientific validation, data requirement, risk, priority, acceptance |
|---|---|---|---|
| TCR mode / TRB | Implemented | adapters and original RFU backend | Bounded public; P0; parity/reference coverage complete. |
| TRA | Partially implemented | schema/adapters only | No TRA RFU reference validation; high; P1; receptor-specific reference required. |
| Paired chain | Blocked | records can coexist; no joint assignment | No method/reference; high; P2; define and validate joint model. |
| BCR data model | Implemented | IG chain normalization/canonical rows | Not a functional model; medium; P2; schema tests only. |
| BCR functional reference, heavy/light, isotype, SHM, family support | Blocked | design gate only | Dedicated reference/public/heldout data required; high; P2; all BCR gates pass. |

## Software and release

| Capability | Status | Public API / implementation / tests / docs | Scientific validation, data requirement, risk, priority, acceptance |
|---|---|---|---|
| Package metadata | Implemented this phase | `pyproject.toml`; metadata tests | Version remains 0.1.0 pending policy; medium; P0; no placeholders, artifact parity. |
| API reference/tutorials | Partially implemented | API freeze and focused docs/examples | Longitudinal/transfer tutorial pending; medium; P0; clean-user walkthrough. |
| CI / Python 3.10–3.12 | Implemented | GitHub workflow | Remote execution must pass; medium; P0. |
| Small public test data | Partially implemented | synthetic fixtures/examples | Redistributable real fixture undecided; low; P1. |
| Versioning | Blocked | changelog/history docs | Existing tag/package inconsistency requires maintainer decision; high; P0. |
| Open-source license | Implemented | MIT `LICENSE` | Low; P0; included in wheel/sdist. |
| DOI archive readiness | Blocked | submission gate | Needs approved versioned release; medium; P1. |
| Reproducibility bundle | Partially implemented | manifests/workflows/plans | Real public manifests/source tables absent; high; P0. |

## Month 1 conclusion

The software foundation now extends beyond orchestration, but the central
transferability claim remains unvalidated on independent real cohorts. The
highest-risk gaps are original-reference compatibility across chain modes,
prespecified longitudinal validation, two independent TCR cohorts including one
held out, and complete bounded scaling/parity tables. BCR remains out of scope.
