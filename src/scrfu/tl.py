from __future__ import annotations

from pathlib import Path
from typing import Any, Optional, Sequence, Union

import pandas as pd

from .tl_rfu_repo import call_rfu_repo

PathLike = Union[str, Path]


def call_rfu(
    adata: Any,
    *,
    backend: str = "rfu_repo",
    rfu_dir: Optional[PathLike] = None,
    chain: str = "TRB",
    airr_key: str = "airr",
    prefer_productive: bool = True,
    wrapper_r_path: PathLike = "r/run_rfu_repo.R",
    rscript_bin: str = "Rscript",
    extra_r_args: Optional[Sequence[str]] = None,
    workdir: Optional[PathLike] = None,
    out_key: str = "rfu",
) -> pd.DataFrame:
    """
    Unified RFU calling entrypoint.

    Parameters
    ----------
    backend
        Currently supported:
          - "rfu_repo": call upstream RFU repo via r/run_rfu_repo.R (requires rfu_dir)
    rfu_dir
        Required when backend="rfu_repo". Path to upstream RFU checkout.
    """
    backend = backend.lower().strip()

    if backend == "rfu_repo":
        if rfu_dir is None:
            raise ValueError("backend='rfu_repo' requires rfu_dir (path to upstream RFU repo checkout).")
        return call_rfu_repo(
            adata,
            rfu_dir=rfu_dir,
            chain=chain,
            airr_key=airr_key,
            prefer_productive=prefer_productive,
            wrapper_r_path=wrapper_r_path,
            rscript_bin=rscript_bin,
            extra_r_args=extra_r_args,
            workdir=workdir,
            out_key=out_key,
        )

    raise ValueError(f"Unknown backend: {backend}. Supported: 'rfu_repo'")


__all__ = ["call_rfu", "call_rfu_repo"]
