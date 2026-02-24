#!/usr/bin/env Rscript

suppressWarnings(suppressMessages({
  library(optparse)
  library(data.table)
}))

option_list <- list(
  make_option(c("--in"), type="character", help="Input TSV with columns: cell_id, cdr3aa, trbv"),
  make_option(c("--out"), type="character", help="Output TSV"),
  make_option(c("--rfu_dir"), type="character", help="Path to upstream RFU repo"),
  make_option(c("--workdir"), type="character", default=NA, help="Scratch directory")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Access using [[ ]] because 'in' is reserved in R
input_path <- opt[["in"]]
output_path <- opt[["out"]]
rfu_dir <- opt[["rfu_dir"]]
workdir <- opt[["workdir"]]

if (is.null(input_path) || is.null(output_path) || is.null(rfu_dir)) {
  stop("Missing required args: --in, --out, --rfu_dir")
}

rfu_dir <- normalizePath(rfu_dir)

rfu_r <- file.path(rfu_dir, "RFU.R")
trimer_rdata <- file.path(rfu_dir, "trimerMDSfit_small.Rdata")
km5000_rdata <- file.path(rfu_dir, "km5000noMax.Rdata")

for (f in c(rfu_r, trimer_rdata, km5000_rdata)) {
  if (!file.exists(f)) stop(paste("Missing required file in rfu_dir:", f))
}

# Workdir default
output_path <- normalizePath(output_path, mustWork = FALSE)
out_dir <- dirname(output_path)

if (is.na(workdir) || is.null(workdir) || nchar(workdir) == 0) {
  workdir <- file.path(out_dir, "work")
}
dir.create(workdir, recursive = TRUE, showWarnings = FALSE)

# Load upstream RFU
source(rfu_r)
load(trimer_rdata)
load(km5000_rdata)

# Read input TSV
inp <- fread(input_path, sep="\t", header=TRUE)

required_cols <- c("cell_id", "cdr3aa", "trbv")
missing <- setdiff(required_cols, names(inp))
if (length(missing) > 0) {
  stop(paste("Input TSV missing columns:", paste(missing, collapse=", ")))
}

# Build GIANA-like TSV
giana <- data.table(
  CDR3 = inp$cdr3aa,
  Vgene = inp$trbv,
  Freq = 1L,
  RANK = 1L,
  Info = inp$cell_id
)

tmp_in <- file.path(workdir, "giana_input.tsv")
fwrite(giana, tmp_in, sep="\t")

# Call RFU
res <- AssignRFUs(tmp_in, CL=km5000noMax)

if (!is.list(res) || is.null(res$TCR) || is.null(res$COR)) {
  stop(paste(
    "Unexpected AssignRFUs output structure.",
    "Expected list with $TCR and $COR.",
    "Got names:", paste(names(res), collapse=", "),
    sep="\n"
  ))
}

tcr <- res$TCR
cor <- res$COR

dt_map <- data.table(
  cdr3aa = names(tcr),
  rfu_id = as.integer(tcr),
  rfu_score = as.numeric(cor[names(tcr)])
)

out <- merge(
  inp,
  dt_map,
  by = "cdr3aa",
  all.x = TRUE,
  sort = FALSE
)

out[, rfu_label := ifelse(is.na(rfu_id), NA_character_, paste0("RFU", rfu_id))]
out <- out[, .(cell_id, cdr3aa, trbv, rfu_label, rfu_score, rfu_id)]

fwrite(out, output_path, sep="\t")
cat("Wrote:", output_path, "\n")
