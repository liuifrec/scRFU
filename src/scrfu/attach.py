from __future__ import annotations

from datetime import datetime, timezone
from typing import Any, Dict, Optional

import pandas as pd

from ._version import __version__


def attach_rfu_results(
    adata: Any,
    features: pd.DataFrame,
    rfu_df: pd.DataFrame,
    *,
    label_col: str = "rfu_label",
    score_col: str = "rfu_score",
    cdr3_col_out: str = "trb_cdr3aa",
    trbv_col_out: str = "trbv",
    uns_key: str = "scrfu",
    provenance: Optional[Dict[str, object]] = None,
) -> None:
    """
    Attach RFU results + extracted TRB features into adata.obs and provenance into adata.uns.

    Requires:
      features: columns cell_id,cdr3aa,trbv
      rfu_df: columns cell_id,rfu_label,rfu_score  (names configurable via args)
    """
    if not hasattr(adata, "obs") or not hasattr(adata, "uns"):
        raise TypeError("adata must look like AnnData (needs .obs and .uns).")

    for col in ["cell_id", "cdr3aa", "trbv"]:
        if col not in features.columns:
            raise ValueError(f"features missing required column: {col}")

    if "cell_id" not in rfu_df.columns:
        raise ValueError("rfu_df must contain 'cell_id'")

    # Merge on cell_id
    feat = features[["cell_id", "cdr3aa", "trbv"]].copy()
    res = rfu_df.copy()

    merged = pd.merge(feat, res, on="cell_id", how="left")

    # Align to obs_names / index
    if hasattr(adata, "obs_names"):
        idx = pd.Index(map(str, list(adata.obs_names)))
    else:
        idx = pd.Index(map(str, list(adata.obs.index)))

    merged = merged.set_index("cell_id").reindex(idx)

    adata.obs[cdr3_col_out] = merged["cdr3aa"].astype("string")
    adata.obs[trbv_col_out] = merged["trbv"].astype("string")

    # Try coerce rfu columns if present
    if label_col in merged.columns:
        adata.obs[label_col] = merged[label_col].astype("string")
    elif "rfu_label" in merged.columns:
        adata.obs[label_col] = merged["rfu_label"].astype("string")
    else:
        # create empty if absent
        adata.obs[label_col] = pd.Series([pd.NA] * len(idx), index=idx, dtype="string")

    if score_col in merged.columns:
        adata.obs[score_col] = pd.to_numeric(merged[score_col], errors="coerce")
    elif "rfu_score" in merged.columns:
        adata.obs[score_col] = pd.to_numeric(merged["rfu_score"], errors="coerce")
    else:
        adata.obs[score_col] = pd.Series([pd.NA] * len(idx), index=idx, dtype="float")

    # Provenance
    prov = {
        "package_version": __version__,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
    }
    if provenance:
        prov.update(provenance)

    adata.uns[uns_key] = prov