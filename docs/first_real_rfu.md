# First real RFU run

Install the software candidate as described in the repository README.
This walkthrough uses small synthetic receptor inputs with **real upstream
RFU assignments**. It is an engineering check, not cohort or radiation-exposure
evidence. BCR outputs from the broader tutorial remain experimental features.

The primary [RFU codebase](https://github.com/s175573/RFU) and
[upstream instructions at the verified revision](https://github.com/s175573/RFU/blob/ff7a1ea6aca28444a7ff8ec7ef09745e7341f99d/README.md)
require `RFU.R`, `trimerMDSfit_small.Rdata` (creates `fit`) and
`km5000noMax.Rdata` (creates `km5000noMax`) in one external directory.
The upstream [paper's code availability statement](https://www.nature.com/articles/s41467-024-52522-z)
identifies this repository. Follow its terms when obtaining these files.

With Git and R installed, obtain a separate checkout and inspect the selected
revision; do not reuse a modified checkout without recording its differences:

```bash
git clone https://github.com/s175573/RFU.git /path/to/external/RFU
git -C /path/to/external/RFU checkout ff7a1ea6aca28444a7ff8ec7ef09745e7341f99d
Rscript -e 'install.packages(c("optparse", "data.table"), repos="https://cloud.r-project.org")'
```

The standard assignment path in this revision uses base R/statistics; scRFU's
wrapper additionally requires `optparse` and `data.table`. Other upstream
analysis functions may have further dependencies and are not invoked here.
The local check used R 4.1.2, optparse 1.7.5 and data.table 1.18.2.1.

```bash
export RFU_DIR=/path/to/external/RFU
scrfu doctor --verbose
python -m scrfu.tutorial --backend rfu_repo --rfu-dir "$RFU_DIR" \
  --outdir /tmp/scrfu-first-real
```

In PowerShell set `$env:RFU_DIR = 'C:\path\to\RFU'` and use a Windows output
path. `Rscript` must be on PATH; the table API and `call-rfu`/`doctor` commands
also accept an explicit `rscript_bin` / `--rscript-bin`. The packaged wrapper is
selected independently of the working directory. A custom wrapper can be passed
with `wrapper_r_path` / `--wrapper-r-path`; doctor accepts the same override.

Expected files include `assigned_receptors.tsv` (12 rows),
`rfu_pseudobulk_proportions.tsv`, `longitudinal_similarity.tsv`, and
`run_manifest.json` with `backend: "rfu_repo"`, wrapper and external-file hashes.
The nearest RFU label remains present even below threshold 0.6; `pass_thr`
distinguishes threshold-qualified assignments. An RFU label alone is not an
antigen-specificity claim. Exact CDR3 duplicates map back to every receptor row;
the standard upstream embedding does not use V-gene differences for assignment.

A separate Linux/macOS parity command compares scRFU directly with untouched
upstream `AssignRFUs(..., CL=km5000noMax, THR=0.6)` on a reordered 12-row fixture:

```bash
python /absolute/path/to/reviewed/scRFU/examples/original_rfu_parity.py \
  --rfu-dir "$RFU_DIR" --outdir /tmp/scrfu-engineering-parity \
  --threshold 0.6 --score-tolerance 1e-12
```

This script imports the installed package and does not inject repository source.
It uses the POSIX `resource` module, so on Windows run this parity script in WSL;
the packaged tutorial is platform-independent Python; actual Windows R execution requires separate verification. Expect `summary.json` with
`passed: true`, 10 eligible rows, 5 distinct eligible CDR3s, 2 ineligible rows,
2 below-threshold rows and zero label, threshold and reconstruction mismatches.
The verified local check had maximum absolute score error 0.0 at tolerance 1e-12.
Use a fresh output directory and retain the manifest with model hashes.
