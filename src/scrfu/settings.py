from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class RFUSettings:
    chain: str = "TRB"
    airr_key: str = "airr"
    out_key: str = "rfu"
    label_col: str = "rfu_label"
    score_col: str = "rfu_score"
    cdr3_col: str = "trb_cdr3aa"
    trbv_col: str = "trbv"