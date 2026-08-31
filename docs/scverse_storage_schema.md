# Native AnnData and MuData storage schema

The stable native storage schema is version `1.0`. For AnnData, source receptor
chains remain unchanged in `adata.obsm["airr"]`; aligned RFU records are stored
in `adata.obsm["scrfu"]`. For MuData, both live in the selected AIRR modality,
normally `mdata.mod["airr"]`.

Each observation contains a variable-length chain list. The AIRR and scRFU
lists must have exactly the same length and order. Every scRFU chain record has:

- `eligible`, `rfu_id`, `rfu_label`, `rfu_score`, and `threshold_pass`;
- `assignment_status` and `unique_sequence_id`;
- `locus`, zero-based `chain_index`, and portable `reference_identity`.

TRA, BCR, nonproductive, missing-sequence and malformed chains remain explicit
non-eligible records. Multiple TRBs are never collapsed. `rfu_id` is the
nearest assignment; `threshold_pass` states whether that assignment satisfies
the configured threshold.

Portable provenance is stored in the AIRR modality's `uns["scrfu"]`. It records
schema/software versions, frozen-reference hashes/identity, mode, threshold,
chain rules, adapter keys, and any explicit summary policy. Absolute runtime
paths are deliberately excluded.

Cell-level `obs` values are opt-in. The default `ambiguity_aware` policy emits
no RFU when eligible chains disagree. Other explicit policies are
`highest_score`, `highest_threshold_score`, `first_eligible`, and
`primary_vdj`. The last requires Scirpy chain indices.

`validate_scrfu_schema` checks alignment, required fields, schema version and
reference identity. `concat_scrfu` validates that every input uses the same
schema and frozen reference before delegating to `anndata.concat`; incompatible
references raise an error. Ordinary `anndata.concat` has no scRFU-specific
reference guard and therefore should not be used to assert reference
compatibility.

H5AD/H5MU round trips and observation subsetting are tested. Awkward support is
currently marked experimental by AnnData, so scRFU versions its own schema and
will document migrations if upstream serialization changes. Backward
compatibility is maintained within a released schema where practical; a future
incompatible schema will require an explicit version change. Stable public API
removals follow deprecation warnings across a release cycle. Experimental BCR
APIs have a weaker compatibility promise and do not share this RFU schema.

Native in-memory execution reads only `obs_names` and the named AIRR `obsm`
entry; it does not inspect `X`, `raw`, layers or `var`. Extremely large H5AD
files should continue to use the targeted receptor-only reader. Full backed
native AIRR mutation is not claimed.

