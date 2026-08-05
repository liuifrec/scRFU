#!/usr/bin/env Rscript

suppressWarnings(suppressMessages({
  library(optparse)
  library(data.table)
}))

option_list <- list(
  make_option(c("--in"), type="character", help="Input TSV with columns: unique_sequence_id, cdr3aa, trbv"),
  make_option(c("--out"), type="character", help="Output TSV"),
  make_option(c("--rfu_dir"), type="character", help="Path to upstream RFU repo"),
  make_option(c("--workdir"), type="character", default=NA, help="Scratch directory"),
  make_option(c("--mode"), type="character", default="standard", help="Backend mode: standard or map_aware"),
  make_option(c("--threshold"), type="double", default=0.6, help="RFU correlation threshold")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Access using [[ ]] because 'in' is reserved in R
input_path <- opt[["in"]]
output_path <- opt[["out"]]
rfu_dir <- opt[["rfu_dir"]]
workdir <- opt[["workdir"]]
mode <- opt[["mode"]]
threshold <- opt[["threshold"]]

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

required_cols <- c("unique_sequence_id", "cdr3aa", "trbv")
missing <- setdiff(required_cols, names(inp))
if (length(missing) > 0) {
  stop(paste("Input TSV missing columns:", paste(missing, collapse=", ")))
}
if (anyDuplicated(inp$unique_sequence_id)) {
  stop("Input TSV contains duplicate unique_sequence_id values")
}
if (any(!grepl("^C", inp$cdr3aa))) {
  stop("Input TSV contains an upstream-ineligible CDR3 that does not start with C")
}

# Build GIANA-like TSV
giana <- data.table(
  CDR3 = inp$cdr3aa,
  Vgene = inp$trbv,
  Freq = 1L,
  RANK = 1L,
  Info = inp$unique_sequence_id
)

tmp_in <- file.path(workdir, "giana_input.tsv")
# EncodeRepertoire() always reads header=FALSE. Writing column names would add a
# synthetic "CDR3" row that itself passes the upstream ^C filter.
fwrite(giana, tmp_in, sep="\t", col.names=FALSE)

# Call the explicitly selected RFU function. Map-aware MAP identifiers are not
# used for reconstruction because the optional extension maps raw metadata rows
# positionally after EncodeRepertoire may have filtered them.
if (identical(mode, "standard")) {
  if (!exists("AssignRFUs", mode="function")) {
    stop("Standard mode requires AssignRFUs(), but it was not found in RFU.R")
  }
  res <- AssignRFUs(tmp_in, CL=km5000noMax, THR=threshold)
} else if (identical(mode, "map_aware")) {
  if (!exists("AssignRFUs_with_map", mode="function")) {
    stop("Map-aware mode requires AssignRFUs_with_map(), but it was not found in RFU.R")
  }
  res <- AssignRFUs_with_map(tmp_in, CL=km5000noMax, THR=threshold)
} else {
  stop(paste("Unknown backend mode:", mode))
}

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
if (length(tcr) != nrow(inp) || length(cor) != nrow(inp)) {
  stop(paste(
    "RFU assignment length does not match submitted query rows.",
    "TCR:", length(tcr), "COR:", length(cor), "input:", nrow(inp)
  ))
}

dt_map <- data.table(
  unique_sequence_id = inp$unique_sequence_id,
  rfu_id = as.integer(tcr),
  rfu_score = as.numeric(cor),
  pass_thr = as.numeric(cor) >= threshold,
  upstream_n_miss = as.integer(res$N)
)

dt_map[, rfu_label := ifelse(is.na(rfu_id), NA_character_, paste0("RFU", rfu_id))]
out <- dt_map[, .(
  unique_sequence_id,
  rfu_id,
  rfu_label,
  rfu_score,
  pass_thr,
  upstream_n_miss
)]

fwrite(out, output_path, sep="\t")
cat("Wrote:", output_path, "\n")
