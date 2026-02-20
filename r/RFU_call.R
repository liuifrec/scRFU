#!/usr/bin/env Rscript

suppressWarnings(suppressMessages({
  library(optparse)
  library(data.table)
}))

option_list <- list(
  make_option(c("--in"), type="character", help="Input TSV with cell_id,cdr3aa,trbv"),
  make_option(c("--atlas"), type="character", help="Atlas path (km5000 centroids)"),
  make_option(c("--out"), type="character", help="Output TSV with cell_id,rfu_label,rfu_score")
)

opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$in) || is.null(opt$atlas) || is.null(opt$out)) {
  stop("Missing required args: --in, --atlas, --out")
}

# NOTE: This is a placeholder scaffold.
# Replace this logic with the real RFU calling implementation.
inp <- fread(opt$in, sep="\t", header=TRUE)

# naive dummy labeling: everything becomes RFU_UNKNOWN with score NA
out <- data.table(
  cell_id = inp$cell_id,
  rfu_label = "RFU_UNKNOWN",
  rfu_score = NA_real_
)

fwrite(out, opt$out, sep="\t")
cat("Wrote:", opt$out, "\n")