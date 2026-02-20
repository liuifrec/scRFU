from __future__ import annotations

from typing import Any, Dict, Optional

import pandas as pd

from ._version import __version__


def _pick_first(df: pd.DataFrame, candidates: list[str]) -> Optional[str]:
    for c in candidates:
        if c in df.columns:
            return c
    return None


def attach_rfu_results(
    adata: Any,
    features: pd.DataFrame,
    rfu_df: pd.DataFrame,
    provenance: Optional[Dict] = None,
    out_key: str = "rfu",
    cdr3_col_out: str = "trb_cdr3aa",
    trbv_col_out: str = "trbv",
    label_col_out: str = "rfu_label",
    score_col_out: str = "rfu_score",
) -> None:
    """
    Attach RFU results back into AnnData-like objects.

    Schema-robust:
      - features must include cell_id and a CDR3aa column (cdr3aa / trb_cdr3aa / cdr3_aa / cdr3)
      - features must include a TRBV column (trbv / v_call / trbv_call / v_gene / vgene)

    Writes:
      - adata.obs[cdr3_col_out], adata.obs[trbv_col_out]
      - adata.obs[label_col_out], adata.obs[score_col_out]
      - adata.uns["scrfu"]["package_version"]
      - adata.uns["scrfu"][out_key] provenance
    """
    if provenance is None:
        provenance = {}

    # Always ensure scrfu uns root exists + package_version (back-compat)
    adata.uns.setdefault("scrfu", {})
    adata.uns["scrfu"]["package_version"] = __version__
    adata.uns["scrfu"].setdefault(out_key, {})

    # Handle empty features (still create obs columns + record provenance)
    if features is None or len(features) == 0:
        if hasattr(adata, "obs"):
            adata.obs[cdr3_col_out] = pd.Series(index=adata.obs.index, dtype="string")
            adata.obs[trbv_col_out] = pd.Series(index=adata.obs.index, dtype="string")
            adata.obs[label_col_out] = pd.Series(index=adata.obs.index, dtype="string")
            adata.obs[score_col_out] = pd.Series(index=adata.obs.index, dtype="float64")

        # Store provenance under the out_key bucket
        adata.uns["scrfu"][out_key].update(provenance)
        return

    if "cell_id" not in features.columns:
        raise ValueError("features must contain column 'cell_id'")

    cdr3_col = _pick_first(features, ["cdr3aa", "trb_cdr3aa", "cdr3_aa", "cdr3"])
    trbv_col = _pick_first(features, ["trbv", "v_call", "trbv_call", "v_gene", "vgene"])

    if cdr3_col is None:
        raise ValueError(f"features is missing a CDR3aa column. Have: {list(features.columns)}")
    if trbv_col is None:
        raise ValueError(f"features is missing a TRBV column. Have: {list(features.columns)}")

    feat = features.copy()
    feat = feat.rename(columns={cdr3_col: "cdr3aa", trbv_col: "trbv"})

    # rfu_df must contain cell_id + label/score (tolerant)
    if rfu_df is None:
        rfu_df = pd.DataFrame(columns=["cell_id", "rfu_label", "rfu_score"])

    if "cell_id" not in rfu_df.columns:
        raise ValueError("rfu_df must contain column 'cell_id'")

    label_col = _pick_first(rfu_df, ["rfu_label", "label", "RFU", "rfu"])
    score_col = _pick_first(rfu_df, ["rfu_score", "score", "COR", "cor", "distance", "dist", "prob"])

    rf = rfu_df.copy()
    if label_col is None:
        rf["rfu_label"] = pd.NA
    else:
        rf = rf.rename(columns={label_col: "rfu_label"})

    if score_col is None:
        rf["rfu_score"] = pd.NA
    else:
        rf = rf.rename(columns={score_col: "rfu_score"})

    merged = feat.merge(rf[["cell_id", "rfu_label", "rfu_score"]], on="cell_id", how="left")

    # Attach by aligning to adata.obs index (assumes obs index are cell IDs)
    by_cell = merged.set_index("cell_id").reindex(adata.obs.index)

    adata.obs[cdr3_col_out] = by_cell["cdr3aa"].astype("string")
    adata.obs[trbv_col_out] = by_cell["trbv"].astype("string")
    adata.obs[label_col_out] = by_cell["rfu_label"].astype("string")
    adata.obs[score_col_out] = pd.to_numeric(by_cell["rfu_score"], errors="coerce")

    # Store provenance in the out_key bucket
    adata.uns["scrfu"][out_key].update(provenance)
