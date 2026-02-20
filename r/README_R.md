# RFU R script interface

`scrfu` calls `Rscript r/RFU_call.R` with the following arguments:

- `--in <input.tsv>` : TSV with columns `cell_id`, `cdr3aa`, `trbv`
- `--atlas <atlas.tsv.gz>` : km5000 atlas / centroids
- `--out <output.tsv>` : TSV output with columns `cell_id`, `rfu_label`, `rfu_score`

You can extend the script to emit additional columns; Python will pass them through.