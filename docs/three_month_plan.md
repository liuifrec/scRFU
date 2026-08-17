# Three-month methods-paper plan

## Month 1 — method and API freeze

Completion criteria:

- capability-by-capability gap analysis and novelty boundary reviewed;
- stable, compatibility, experimental, internal, and deferred API symbols
  recorded;
- original-RFU parity acceptance criteria and artifact requirements frozen;
- AIRR/scirpy/AnnData/optional-MuData contracts documented;
- deterministic parallel chunk implementation passes fake-runner parity;
- reference-coverage terminology and formulas tested;
- repeated-measures validator, longitudinal representations, similarity,
  retrieval, trajectory classification, and donor-block resampling pass
  hand-calculable synthetic tests;
- frozen-reference transfer, explicit harmonization, held-out manifest,
  robustness metrics, and comparator interfaces are synthetically tested;
- BCR remains behind a documented go/no-go gate;
- benchmark dataset roles, figures, source tables, and submission gates are
  frozen without claiming missing validation.

Month 1 does not require private-cohort execution or complete public cohorts.

## Month 2 — benchmarks and applications

Completion criteria:

- row/score/threshold parity against the frozen original RFU artifacts;
- 1k/10k/25k bounded scaling, memory, cache reuse, chunk-size, and serial/parallel
  measurements;
- deterministic downsampling and order-perturbation analyses;
- repeated-measures methods demonstration on the separately governed deep
  longitudinal cohort, framed without population-level claims;
- at least one independent public TCR cohort and one public paired single-cell
  application analyzed under frozen definitions;
- held-out cohort and evaluation metrics registered before evaluation;
- comparator performance reported through the same sample-level interfaces;
- BCR proceeds only if its gate is passed.

## Month 3 — manuscript and release

Completion criteria:

- figures and complete source tables regenerate from manifests;
- tutorials and a small redistributable public/synthetic test dataset pass from
  a clean installation;
- CI, package metadata, license, API documentation, and supported Python matrix
  pass;
- version/tag/release consistency is resolved and artifacts are archived with a
  DOI;
- manuscript, methods, cover letter, and reproducibility report align with the
  evidence actually produced;
- reviewer-shareable code uses public data or documented user-supplied inputs
  and contains no private cohort data.
