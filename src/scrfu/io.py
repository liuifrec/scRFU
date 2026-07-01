from __future__ import annotations

from pathlib import Path

import anndata as ad
import pandas as pd

from .summary import aggregate_rfu

PathLike = str | Path


def read_h5ad(path: PathLike) -> ad.AnnData:
    return ad.read_h5ad(str(path))


def write_h5ad(adata: ad.AnnData, path: PathLike) -> None:
    adata.write_h5ad(str(path))


def export_rfu_matrix(
    adata: object,
    groupby: str,
    path: PathLike,
    normalize: bool = True,
    sep: str = "\t",
) -> pd.DataFrame:
    """
    Export a group-by-RFU matrix to TSV or CSV.
    """
    matrix = aggregate_rfu(adata, groupby=groupby, normalize=normalize)
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    matrix.to_csv(output_path, sep=sep)
    return matrix
