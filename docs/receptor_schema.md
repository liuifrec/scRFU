# Canonical Receptor Schema

Schema version `1.0` represents one receptor chain or selected record per row.
Required columns are `input_row_id`, `cell_id`, `chain`, `cdr3aa`, `v_call`,
`productive`, `source_adapter`, and `source_row_id`.

`input_row_id` is deterministic and unique within a prepared table. `cell_id`
is a string. Recognized loci are `TRA`, `TRB`, `TRG`, `TRD`, `IGH`, `IGK`, and
`IGL`. CDR3 and V-call text retain source information after safe trimming and
missing-value normalization. V calls are provenance, not RFU deduplication
keys. Source order remains stable.

Standard optional fields are `duplicate_count`, `consensus_count`,
`source_slot`, `j_call`, `d_call`, `c_call`, `junction`, `junction_aa`,
`umi_count`, `read_count`, `clonotype_id`, `receptor_type`, and `light_chain`.
Cell metadata remain separate by default.

`scrfu.pp.validate_receptor_table()` reports counts, identifier uniqueness,
missing identifiers and CDR3s, chains, productive normalization, duplicated
source IDs, schema version, status, and actionable errors. Canonical validity
does not require RFU eligibility: non-TRB and non-C-starting records are valid.
