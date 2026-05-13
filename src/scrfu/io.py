from __future__ import annotations

from pathlib import Path

import anndata as ad

PathLike = str | Path


def read_h5ad(path: PathLike) -> ad.AnnData:
    return ad.read_h5ad(str(path))


def write_h5ad(adata: ad.AnnData, path: PathLike) -> None:
    adata.write_h5ad(str(path))
