# Draft scverse ecosystem submission

This is a non-submitted working draft. `docs/scverse_meta_draft.yaml` follows
registry schema 2.0 and uses the controlled primary category `Adaptive immune
cell receptor` with accepted tags only.

Proposed description of data-structure use: scRFU reads current Scirpy AIRR
Awkward records from AnnData or a MuData AIRR modality, writes positionally
aligned chain-level RFU records to `obsm`, and stores portable schema/reference
provenance in `uns`. Optional cell summaries are explicitly requested and are
written to modality `obs`. RFU execution is independent of expression `X`.

Mandatory criteria not yet met:

- Standard registry: `install.pypi: scrfu` is prospective until the maintainer
  publishes and validates the distribution. The draft must not be submitted
  before that identifier resolves to this project.
- Hosted API docs: the URL is prospective until the maintainer enables a docs
  service and verifies `objects.inv`.
- Maintainer consent: the maintainer must explicitly agree to listing in the
  eventual PR.

The future PR checklist must also state that official RFU code/reference assets
are external and that Python-only interoperability CI does not establish
official RFU execution on macOS or Windows.

