#!/usr/bin/env Rscript

suppressWarnings(suppressMessages({
  library(optparse)
  library(data.table)
}))

option_list <- list(
  make_option(c("--in"), type="character", help="Input TSV: cell_id, cdr3aa, trbv"),
  make_option(c("--out"), type="character", help="Output TSV: cell_id, rfu_label, rfu_score"),
  make_option(c("--rfu_dir"), type="character", help="Path to upstream RFU repo (contains RFU.R and Rdata)")
)

opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$in) || is.null(opt$out) || is.null(opt$rfu_dir)) {
  stop("Missing required args: --in, --out, --rfu_dir")
}

rfu_dir <- normalizePath(opt$rfu_dir)

rfu_r <- file.path(rfu_dir, "RFU.R")
trimer_rdata <- file.path(rfu_dir, "trimerMDSfit_small.Rdata")
km5000_rdata <- file.path(rfu_dir, "km5000noMax.Rdata")

for (f in c(rfu_r, trimer_rdata, km5000_rdata)) {
  if (!file.exists(f)) {
    stop(paste("Missing required file in rfu_dir:", f))
  }
}

# Load upstream RFU code + atlas objects (names are defined upstream)
source(rfu_r)
load(trimer_rdata)
load(km5000_rdata)

# Read input features
inp <- fread(opt$in, sep="\t", header=TRUE)

required_cols <- c("cell_id", "cdr3aa", "trbv")
missing <- setdiff(required_cols, names(inp))
if (length(missing) > 0) {
  stop(paste("Input TSV missing columns:", paste(missing, collapse=", ")))
}

# ----------------------------
# Upstream RFU call strategy
# ----------------------------
# The upstream repo provides functions in RFU.R.
# We intentionally keep this wrapper minimal and let upstream define details.
#
# We attempt a best-effort mapping:
# - Use cdr3aa as the sequence input
# - Optionally use TRBV as gene context if supported by upstream functions
#
# Because upstream function names may evolve, we try a small set of known options.
# If none match, error with a clear message.

assign_ok <- FALSE
out <- NULL

# Option 1: AssignRFUs (as mentioned in upstream usage patterns)
if (exists("AssignRFUs")) {
  # Common pattern: AssignRFUs(seq, trimerMDSfit_small, km5000noMax, ...)
  # We'll try different argument forms robustly.
  try({
    out <- AssignRFUs(inp$cdr3aa)
    assign_ok <- TRUE
  }, silent = TRUE)

  if (!assign_ok) {
    try({
      out <- AssignRFUs(inp$cdr3aa, inp$trbv)
      assign_ok <- TRUE
    }, silent = TRUE)
  }
}

# Option 2: RFUassign / assignRFU (fallback guesses)
if (!assign_ok && exists("RFUassign")) {
  try({
    out <- RFUassign(inp$cdr3aa)
    assign_ok <- TRUE
  }, silent = TRUE)
}
if (!assign_ok && exists("assignRFU")) {
  try({
    out <- assignRFU(inp$cdr3aa)
    assign_ok <- TRUE
  }, silent = TRUE)
}

if (!assign_ok) {
  stop(paste(
    "Could not find a usable RFU assignment function in RFU.R.",
    "Tried: AssignRFUs, RFUassign, assignRFU.",
    "Please check upstream RFU repo API and update r/run_rfu_repo.R accordingly.",
    sep = "\n"
  ))
}

# ----------------------------
# Normalize output
# ----------------------------
# We expect out to be a data.frame-like with per-sequence RFU labels/scores.
# Since upstream output schema may vary, we normalize to:
#   cell_id, rfu_label, rfu_score
#
# Heuristics:
# - if out has column 'RFU' or 'rfu' -> label
# - else if out has column 'label' -> label
# - score: 'score', 'dist', 'distance', 'prob' (first found)
#
# If out is a vector, treat it as label only.

if (is.vector(out) && !is.list(out)) {
  dt <- data.table(cell_id = inp$cell_id, rfu_label = as.character(out), rfu_score = NA_real_)
} else {
  dt_out <- as.data.table(out)

  # Find label column
  label_col <- NULL
  for (c in c("RFU", "rfu", "label", "cluster", "centroid")) {
    if (c %in% names(dt_out)) { label_col <- c; break }
  }
  if (is.null(label_col)) {
    stop(paste("RFU output missing label-like column. Columns:", paste(names(dt_out), collapse=", ")))
  }

  score_col <- NULL
  for (c in c("score", "dist", "distance", "prob", "p", "similarity")) {
    if (c %in% names(dt_out)) { score_col <- c; break }
  }

  if (is.null(score_col)) {
    dt <- data.table(
      cell_id = inp$cell_id,
      rfu_label = as.character(dt_out[[label_col]]),
      rfu_score = NA_real_
    )
  } else {
    dt <- data.table(
      cell_id = inp$cell_id,
      rfu_label = as.character(dt_out[[label_col]]),
      rfu_score = as.numeric(dt_out[[score_col]])
    )
  }
}

fwrite(dt, opt$out, sep="\t")
cat("Wrote:", opt$out, "\n")