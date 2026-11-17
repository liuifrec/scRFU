#!/usr/bin/env Rscript
# Compatibility launcher for source-checkout users; the canonical wrapper is packaged.
file_arg <- grep("^--file=", commandArgs(trailingOnly=FALSE), value=TRUE)
if (length(file_arg) != 1L) stop("Run this launcher with Rscript")
launcher <- normalizePath(sub("^--file=", "", file_arg[[1]]))
source(file.path(dirname(launcher), "..", "src", "scrfu", "resources", "run_rfu_repo.R"))
