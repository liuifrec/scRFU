# scRFU 0.4.0 release notes

These notes describe the prepared 0.4.0 source. They do not indicate that a Git
tag, package upload, hosted release, archive, or DOI has been created.

## Overview

scRFU 0.4.0 provides deterministic, frozen-reference TRB RFU assignment and
receptor-level analysis across table, AIRR/scirpy, AnnData, MuData-like, Cell
Ranger VDJ, and Wells inputs. It preserves the official upstream RFU standard
mode while adding exact-CDR3 query deduplication, stable input-row
reconstruction, reference-coverage diagnostics, and portable provenance.

RFU execution uses a user-supplied external RFU checkout and assets; they are
not bundled. Official RFU/R execution has been validated on the primary Linux
development environment. Python-only package functionality is CI-tested on
Linux, macOS, and Windows with Python 3.10, 3.11, and 3.12 on the Linux matrix.
This does not claim official RFU/R execution support on macOS or Windows.

## Analysis and execution

- Restartable, checksummed chunk execution supports deterministic serial and
  parallel orchestration, exact cache validation, and rapid reuse of completed
  work.
- Selective H5AD and receptor-cache paths avoid expression-matrix loading and
  support bounded-memory receptor processing at Wells-atlas scale.
- Analysis utilities cover conventional repertoire metrics, RFU abundance and
  convergence, pseudobulk, overlap, phenotype coupling, repeated-measures
  design validation, longitudinal similarity/retrieval, frozen-reference
  transfer, and deterministic robustness/downsampling.
- The optional VDJdb framework consumes an explicit, version-labelled local
  reference. It performs exact CDR3 or strict CDR3+V evidence matching and
  annotation-coherence comparisons without bundling or silently downloading
  VDJdb. Matches are external annotation evidence, not proof of antigen
  specificity.
- `scrfu doctor` reports installation and external-backend configuration while
  redacting full paths unless `--verbose` is requested.

## Documentation and distribution

The distribution includes a small synthetic fixture and runnable end-to-end
tutorial. Core installation requires AnnData, NumPy, and pandas. Plotting,
scirpy, and MuData integrations remain optional extras. Sphinx sources document
the stable TCR API separately from experimental BCR preprocessing.

scRFU code is MIT licensed. Upstream RFU code/assets, RFU reference `.Rdata`,
VDJdb, Wells, GEO, OAS, 10x, and other public datasets are separate external
resources and are not redistributed.

## Experimental BCR scope

0.4.0 includes experimental BCR canonicalization, QC, deterministic heavy/light
pairing, and conservative feature extraction. It does **not** include a BCR
functional-unit reference or assignment API. The public-data feasibility gate
remained NO-GO, and TCR RFUs must not be applied to IGH/IGK/IGL receptors.

## Known limitations

- Users must obtain and configure the upstream RFU implementation and reference
  assets separately for real RFU assignment.
- Full external RFU/R execution is validated on the primary Linux environment;
  cross-platform CI covers Python-only functionality.
- VDJdb analyses require a separately obtained, explicitly pinned local
  reference and depend on its annotation coverage.
- Threshold failure indicates low similarity to the frozen RFU reference. It is
  not calibrated confidence and is not by itself an out-of-distribution call.
- BCR functionality remains preprocessing and state-feature extraction only.
