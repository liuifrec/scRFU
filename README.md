# scRFU (scrfu)

RFU calling for 10x single-cell TCR data integrated into scirpy/scverse
AnnData.

## Overview

scRFU extracts per-cell TRB features (CDR3aa + TRBV), runs an RFU R
script (km5000 centroids atlas), and attaches RFU assignments back into
AnnData.

## Installation

Editable install for development:

``` bash
pip install -e ".[dev]"
```

Optional scirpy extras:

``` bash
pip install -e ".[dev,scirpy]"
```

## Quickstart (Python)

``` python
import anndata as ad
import scrfu as rf

adata = ad.read_h5ad("input.h5ad")

rf.tl.call_rfu(
    adata,
    chain="TRB",
    airr_key="airr",
    out_key="rfu",
    rscript_path="r/RFU_call.R",
    atlas_path="data/km5000_centroids.tsv.gz",
)

adata.write_h5ad("out.h5ad")
```

## AnnData Contract

### adata.obs

-   rfu_label
-   rfu_score
-   trb_cdr3aa
-   trbv

### adata.uns\["scrfu"\]

Stores provenance: package version, atlas hash, R script hash,
timestamp.
