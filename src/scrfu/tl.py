from __future__ import annotations

from pathlib import Path
from typing import Any, Optional, Sequence, Union

import pandas as pd

from .attach import attach_rfu_results
from .extract import extract_trb_features
from .rfu import file_sha256, run_rfu_rscript


PathLike = Union[str, Path]


def call_rfu(
    adata: Any,
    *,
    chain: str = "TRB",
    airr_key: str = "airr",
    out_key: str = "rfu",
    rscript_path: PathLike,
    atlas_path: PathLike,
    extra_r_args: Optional[Sequence[str]] = None,
) -> pd.DataFrame:
    """
    End-to-end RFU calling:
      1) Extract TRB features per cell from AIRR-like table in adata.obsm[airr_key]
      2) Call RFU via Rscript
      3) Attach results to adata.obs and provenance to adata.uns["scrfu"]

    Returns the RFU result dataframe (cell_id, rfu_label, rfu_score, ...).
    """
    features = extract_trb_features(adata, airr_key=airr_key, chain=chain)

    # if no TRB found, attach minimal provenance and return empty
    if features.empty:
        attach_rfu_results(
            adata,
            features=features,
            rfu_df=pd.DataFrame({"cell_id": [], "rfu_label": [], "rfu_score": []}),
            provenance={
                "chain": chain,
                "airr_key": airr_key,
                "out_key": out_key,
                "note": "No features extracted; empty TRB set.",
            },
        )
        return pd.DataFrame({"cell_id": [], "rfu_label": [], "rfu_score": []})

    out_path = Path(f"{out_key}.tsv")
    run = run_rfu_rscript(
        features,
        rscript_path=rscript_path,
        atlas_path=atlas_path,
        out_path=out_path,
        extra_args=list(extra_r_args or []),
    )

    provenance = {
        "chain": chain,
        "airr_key": airr_key,
        "out_key": out_key,
        "rscript_path": str(rscript_path),
        "atlas_path": str(atlas_path),
        "rscript_sha256": file_sha256(rscript_path) if Path(rscript_path).exists() else None,
        "atlas_sha256": file_sha256(atlas_path) if Path(atlas_path).exists() else None,
    }

    attach_rfu_results(adata, features=features, rfu_df=run.df, provenance=provenance)
    return run.df