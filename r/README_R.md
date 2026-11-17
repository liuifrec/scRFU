# scRFU R wrapper

The canonical scRFU-owned wrapper is `src/scrfu/resources/run_rfu_repo.R` and is
included in wheels and sdists. Python resolves it through installed package
resources; no current-working-directory setup is required. `r/run_rfu_repo.R`
is a compatibility launcher for a source checkout. The historical `RFU_call.R`
script is not the supported backend.

Supply upstream `RFU.R`, `trimerMDSfit_small.Rdata` and `km5000noMax.Rdata`
externally. The wrapper requires R, optparse and data.table. See
[the real RFU tutorial](../docs/first_real_rfu.md) and
[configuration](../docs/configuration.md) for setup, overrides and diagnostics.
