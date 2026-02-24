from __future__ import annotations

from pathlib import Path
from typing import Any, Optional, Sequence, Union

import pandas as pd

from .attach import attach_rfu_results
from .extract import extract_trb_features
from .backends.rfu_repo import RFURepoBackend

PathLike = Union[str, Path]


def call_rfu_repo(
    adata: Any,
    *,
    rfu_dir: PathLike,
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
    Call upstream RFU (https://github.com/s175573/RFU) through scrfu's wrapper and attach results to AnnData.

    Writes:
      - adata.obs["trb_cdr3aa"], adata.obs["trbv"]
      - adata.obs["rfu_label"], adata.obs["rfu_score"]
      - adata.uns["scrfu"] provenance (merged)

    Returns the result dataframe (at least: cell_id, rfu_label, rfu_score).
    """
    features = extract_trb_features(
        adata,
        airr_key=airr_key,
        chain=chain,
        prefer_productive=prefer_productive,
    )

    backend = RFURepoBackend(
        rfu_dir=rfu_dir,
        wrapper_r_path=wrapper_r_path,
        rscript_bin=rscript_bin,
    )

    if features.empty:
        # Attach empty columns + provenance and return
        attach_rfu_results(
            adata,
            features=features,
            rfu_df=pd.DataFrame({"cell_id": [], "rfu_label": [], "rfu_score": []}),
            provenance={
                **backend.provenance_dict(),
                "chain": chain,
                "airr_key": airr_key,
                "out_key": out_key,
                "note": "No TRB features extracted; RFU not run.",
            },
        )
        return pd.DataFrame({"cell_id": [], "rfu_label": [], "rfu_score": []})

    run = backend.run(
        features,
        extra_args=extra_r_args,
        workdir=workdir,
    )

    # Attach results + provenance
    attach_rfu_results(
        adata,
        features=features,
        rfu_df=run.df,
        provenance={
            **backend.provenance_dict(),
            "chain": chain,
            "airr_key": airr_key,
            "out_key": out_key,
        },
    )
    return run.df
