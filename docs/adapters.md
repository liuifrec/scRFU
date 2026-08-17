# Receptor Adapters

Adapters are lightweight callables returning canonical receptors, separate
cell metadata, QC, adapter name/version, and provenance. Use
`list_receptor_adapters()`, `get_receptor_adapter()`, or `prepare_receptors()`.
Formats are selected explicitly.

- `wells_tcr_ir` (`wells`) reads flattened Wells VDJ 1/2 fields and selects a
  primary chain by duplicate count, consensus count, then source slot.
- `scirpy_airr` reads pandas AIRR data from `adata.obsm["airr"]` or directly.
- `generic_airr_dataframe` (`airr`, `dataframe`) accepts AIRR-like DataFrames.
- `cellranger_vdj` (`cellranger`, `tenx_vdj`) accepts Cell Ranger V(D)J
  contig CSV/TSV files or DataFrames without scanpy/scirpy. It maps barcode,
  chain, amino-acid/nucleotide CDR3, V/D/J/C calls, UMI/read counts, and raw
  clonotype/consensus identifiers when present. Primary selection ranks
  productive, confidence, UMI, read count, then stable source order.

AIRR aliases cover cell, chain, CDR3 amino acid, V gene, productive status,
UMI/read, duplicate, and consensus fields. Generic primary selection ranks
productive first, then available UMI, read, duplicate, consensus counts, then
stable order. A priority is used only when its field exists.

Scirpy is optional. Modern awkward-array AIRR objects should currently be
converted to pandas. Cell Ranger requires barcode, chain, and amino-acid CDR3;
all other source columns are optional, including V gene. Contigs with missing
required values are excluded and counted in adapter QC; a missing required
column raises immediately.
