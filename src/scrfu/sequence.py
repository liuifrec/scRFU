from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import pandas as pd

from .repertoire import analysis_frame

AMINO_ACIDS = tuple("ACDEFGHIKLMNPQRSTVWY")


def _visual_alignment(sequence: str, width: int, mode: str) -> str:
    gap = width - len(sequence)
    if mode == "left":
        return sequence + "-" * gap
    if mode == "right":
        return "-" * gap + sequence
    if mode == "center":
        left = gap // 2
        return "-" * left + sequence + "-" * (gap - left)
    if sequence.startswith("C") and sequence.endswith(("F", "W")) and len(sequence) >= 2:
        return sequence[:-1] + "-" * gap + sequence[-1]
    left = gap // 2
    return "-" * left + sequence + "-" * (gap - left)


def rfu_sequence_matrix(
    data: Any,
    *,
    rfu: Any | None = None,
    sequence_col: str = "cdr3aa",
    rfu_col: str = "rfu_label",
    weighting: str = "unique_sequence",
    cell_col: str = "cell_id",
    clonotype_col: str = "clonotype_id",
    alignment: str = "conserved_ends",
    include_gap: bool = True,
) -> pd.DataFrame:
    """Return a deterministic position-frequency matrix for RFU visualization.

    The alignment choices are display conventions, not formal multiple-sequence
    alignment algorithms. ``conserved_ends`` anchors an initial C and terminal
    F/W when both are present and centers other sequences.
    """
    if weighting not in {"unique_sequence", "cell", "clonotype"}:
        raise ValueError("weighting must be 'unique_sequence', 'cell', or 'clonotype'.")
    if alignment not in {"left", "right", "center", "conserved_ends"}:
        raise ValueError("alignment must be 'left', 'right', 'center', or 'conserved_ends'.")
    frame = analysis_frame(data)
    if sequence_col not in frame:
        raise ValueError(f"Sequence matrix input is missing column {sequence_col!r}.")
    if rfu is not None:
        if rfu_col not in frame:
            raise ValueError(f"RFU selection requires column {rfu_col!r}.")
        frame = frame.loc[frame[rfu_col].eq(rfu)].copy()
    frame = frame.loc[frame[sequence_col].notna()].copy()
    frame[sequence_col] = frame[sequence_col].astype(str).str.strip().str.upper()
    frame = frame.loc[frame[sequence_col].ne("")]
    invalid = frame.loc[
        ~frame[sequence_col].map(lambda value: set(value).issubset(AMINO_ACIDS)), sequence_col
    ].drop_duplicates()
    if len(invalid):
        raise ValueError(f"Invalid amino-acid sequences: {invalid.tolist()[:5]}")
    if weighting == "unique_sequence":
        weighted = frame.drop_duplicates(sequence_col).assign(_weight=1.0)
    elif weighting == "cell":
        if cell_col not in frame:
            raise ValueError(f"Cell weighting requires column {cell_col!r}.")
        weighted = frame.drop_duplicates([cell_col, sequence_col]).assign(_weight=1.0)
    else:
        if clonotype_col not in frame or frame[clonotype_col].isna().any():
            raise ValueError(
                f"Clonotype weighting requires complete column {clonotype_col!r}; no sequence fallback is implicit."
            )
        weighted = frame.drop_duplicates([clonotype_col, sequence_col]).copy()
        sequences_per_clonotype = weighted.groupby(clonotype_col, observed=True)[
            sequence_col
        ].transform("nunique")
        weighted["_weight"] = 1.0 / sequences_per_clonotype
    if weighted.empty:
        raise ValueError("No valid sequences are available for the RFU sequence matrix.")
    width = int(weighted[sequence_col].str.len().max())
    alphabet: Sequence[str] = (*AMINO_ACIDS, "-") if include_gap else AMINO_ACIDS
    counts = pd.DataFrame(0.0, index=pd.RangeIndex(1, width + 1, name="position"), columns=alphabet)
    for sequence, weight in weighted.groupby(sequence_col, sort=True)["_weight"].sum().items():
        aligned = _visual_alignment(sequence, width, alignment)
        for position, residue in enumerate(aligned, start=1):
            if residue in counts.columns:
                counts.loc[position, residue] += float(weight)
    denominator = counts.sum(axis=1).replace(0, pd.NA)
    return counts.div(denominator, axis=0).fillna(0.0)


__all__ = ["AMINO_ACIDS", "rfu_sequence_matrix"]
