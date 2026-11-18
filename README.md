# scRFU

scRFU connects frozen-reference TCR repertoire functional units (RFUs) to
single-cell annotation, longitudinal repertoire comparison, and cellular and
regulatory interpretation. It retains receptor rows, threshold status and
provenance when linking assignments to AnnData, Scirpy or MuData. Sample-level
summaries support repeated-measures comparisons without imputing missing visits.
User-supplied RFU-QTL, eQTL, caQTL and GWAS tables provide descriptive evidence
integration; shared variants do not establish colocalization or causality.

The software version is **0.4.1**. Its version-tag installation route is shown
below; the tag is published with the verified GitHub Release. During release
preparation, install the local source with `python -m pip install /path/to/scRFU`.
With Git installed, the version-tag command is:

```bash
python -m venv /tmp/scrfu-env
source /tmp/scrfu-env/bin/activate
python -m pip install "git+https://github.com/liuifrec/scRFU.git@v0.4.1"
scrfu --version
scrfu doctor
python -m scrfu.tutorial --outdir /tmp/scrfu-synthetic
```

On Windows use `py -m venv .venv` and `.venv\Scripts\Activate.ps1`, with Windows output paths. The tutorial works outside the checkout and requires no R or
real data. Its 12 receptor rows and RFU labels are **synthetic mock assignments**;
`run_manifest.json` records this explicitly. See the [executable quickstart](docs/tutorial.md)
for expected outputs and [first real RFU run](docs/first_real_rfu.md) for R setup.

The core package requires Python 3.10+ and pandas 2 or 3; see the
[tested dependency matrix and software status](docs/current_status.md).
Matplotlib, Scirpy/Scanpy and MuData remain optional. Install them from the same
source with `python -m pip install "scrfu[all] @ git+https://github.com/liuifrec/scRFU.git@v0.4.1"`
when needed. For a local checkout, `python -m pip install /path/to/scRFU` also works.

RFU assignment uses a separately supplied [upstream RFU checkout](https://github.com/s175573/RFU).
The scRFU wrapper is packaged; upstream `RFU.R` and its model files remain
external. Set `RFU_DIR` or pass `rfu_dir` / `--rfu-dir`. `scrfu doctor --verbose`
reports the selected wrapper, hashes, R executable and external assets.

- [Single-cell workflow](docs/scverse_tutorial.md)
- [Configuration and wrapper overrides](docs/configuration.md)
- [Regulatory evidence and missing-value contract](docs/regulatory_triangulation.md)
- [Current status and historical records](docs/current_status.md)
- [Release notes and publication procedure](docs/release_notes_0.4.1.md)

Standard RFUs here are TRB-specific. BCR preprocessing is experimental and does
not assign TCR RFUs to BCRs. Reference coverage, negative results and unresolved
reference identities must remain visible in interpretation.

Cite scRFU using [CITATION.cff](CITATION.cff), and cite the upstream method:
Hu J, Pan M, Reid B, Tworoger S, Li B. *Quantifiable blood TCR repertoire
components associate with immune aging*. Nature Communications 15, 8171 (2024).
[doi:10.1038/s41467-024-52522-z](https://doi.org/10.1038/s41467-024-52522-z).
Package-registry publication, a software archive DOI and hosted documentation
are separate publication steps; no unverified endpoint is claimed here. The MIT license
covers scRFU-owned code and synthetic fixtures; see the
[external-resource boundary](docs/release_licensing_boundary.md).
