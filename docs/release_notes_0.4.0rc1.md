# scRFU 0.4.0rc1 release-candidate notes

This is a development release candidate. No tag, package upload, archive, DOI,
or publication has been created.

## TCR RFU framework

- Standardized AIRR/scirpy, AnnData, MuData-like, Cell Ranger, Wells, and
  generic receptor inputs while preserving canonical TRB RFU definitions.
- Added exact-CDR3 query deduplication with stable per-row reconstruction,
  restartable checksummed chunks, deterministic serial/process/thread
  execution, and explicit nearest-versus-threshold-qualified semantics.
- Added reference-coverage diagnostics, conventional repertoire metrics, RFU
  metrics, pseudobulk, overlap, phenotype coupling, longitudinal utilities,
  frozen-reference transfer, comparator utilities, and deterministic robustness
  benchmarks.
- Added long-format, ambiguity-aware VDJdb evidence matching. RFU identity
  remains CDR3-based; strict CDR3+V queries use distinct chain/CDR3/V identities.
  Antigen-label coherence is external annotation evidence, not antigen
  specificity.

## Scale and reproducibility

- Validated a full 610,429-cell Wells receptor-only run: 303,088 productive
  primary-TRB rows, 192,675 unique CDR3 queries, 4,996 nearest RFUs, and 77.18%
  threshold coverage. Runtime outputs and public references remain external.
- Added a generic completed-run validator, immutable external-evidence index,
  public API snapshot, advisory versioned performance baseline, installation
  doctor, synthetic end-to-end fixture/tutorial, and failure-injection tests.
- Added Linux Python 3.10/3.11/3.12 CI intent, Python-only smoke jobs for Linux,
  macOS, and Windows, isolated wheel/sdist checks, and warning-clean Sphinx API
  documentation.

## Experimental BCR utilities

- Added experimental BCR canonicalization, chain/isotype normalization,
  deterministic heavy/light selection, pair reconstruction, QC, conservative
  state features, and a one-row-per-cell feature matrix with explicit
  missingness.
- The public-data feasibility gate is **NO-GO** for a frozen BCR receptor-state
  reference. Exact heavy and paired baselines had 0–0.236% frozen cross-dataset
  coverage; acquired tables lacked SHM and inferred clonal-family fields. No
  BCR functional-unit reference or assignment API is exposed.

## Distribution boundaries

scRFU redistributes only its own code and synthetic fixtures. The upstream RFU
implementation/assets, VDJdb, Wells atlas, GEO/OAS/10x data, and derived
real-data outputs must be obtained and licensed separately by users. Their
hashes and provenance can be recorded without bundling their contents.

