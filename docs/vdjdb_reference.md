# Versioned VDJdb Reference Contract

scRFU uses a user-supplied local table as external annotation evidence. It does
not bundle VDJdb, download a current release, or resolve a changing `latest`
resource. Every load requires a non-empty release label; reproducible analyses
should also pass the expected SHA256 digest.

```python
reference = scrfu.tl.load_vdjdb_reference(
    "vdjdb-release.tsv.gz",
    release_label="EXPLICIT_RELEASE_LABEL",
    expected_sha256="...",
)
```

The input may be a TSV, CSV, gzip-compressed file, or pandas DataFrame. File
delimiter inference recognizes `.tsv` and `.csv` in compressed filenames; `sep`
can override it. `column_mappings` maps canonical names to exact source column
names when aliases are insufficient. Strict validation is the default.

## Canonical reference schema

Only `cdr3aa` is required. Canonical optional fields are `v_call`, `j_call`,
`chain`, `paired_receptor_id`, `paired_cdr3aa`, `paired_v_call`, `epitope`,
`antigen_gene`, `antigen_species`, `mhc`, `mhc_class`, `evidence_score`,
`publication_id`, and `database_row_id`. `reference_row_id` is a deterministic
zero-padded identifier for the input row. Duplicate database records are
reported but are not collapsed, because records can carry distinct evidence.

TRA, TRB, TRG, and TRD are recognized TCR chains. Other chains can remain in
the canonical reference for provenance, but RFU antigen-coherence analysis in
this release is TCR-focused. A missing source chain remains missing; it is not
inferred from the CDR3.

CDR3 normalization consists only of string conversion, whitespace trimming,
uppercasing, and missing-value handling. Invalid amino-acid characters are an
error in strict validation. V genes support `exact` or `strip_allele`; the
latter removes only a terminal allele suffix such as `*01` and never collapses
distinct genes.

## Provenance

`VDJdbReference.provenance` records the explicit release, optional source URL,
runtime filename/path, byte size, SHA256, UTC load time, row and unique-CDR3
counts, available chains and antigen fields, resolved column mappings, scRFU
version, and reference schema version. Runtime manifests may contain the local
path. Static package files do not.

A DataFrame reference receives a SHA256 over a deterministic CSV serialization
and the filename `<dataframe>`. For durable manuscript provenance, a pinned
local file and its expected digest are preferred.

