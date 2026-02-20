from __future__ import annotations

import hashlib
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence, Union

import pandas as pd


PathLike = Union[str, Path]


@dataclass(frozen=True)
class RFURunResult:
    df: pd.DataFrame
    stdout: str
    stderr: str
    returncode: int


def file_sha256(path: PathLike) -> str:
    p = Path(path)
    h = hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def run_rfu_rscript(
    features: pd.DataFrame,
    *,
    rscript_path: PathLike,
    atlas_path: PathLike,
    out_path: Optional[PathLike] = None,
    extra_args: Optional[Sequence[str]] = None,
) -> RFURunResult:
    """
    Call an external R script to perform RFU assignment.

    The R script is expected to accept:
      --in <input.tsv>
      --atlas <atlas.tsv.gz or similar>
      --out <output.tsv>

    Input TSV columns required:
      cell_id, cdr3aa, trbv

    Output TSV expected columns:
      cell_id, rfu_label, rfu_score
    """
    extra_args = list(extra_args or [])

    required = {"cell_id", "cdr3aa", "trbv"}
    if not required.issubset(set(features.columns)):
        raise ValueError(f"features must contain columns {required}; got {list(features.columns)}")

    rscript_path = Path(rscript_path)
    atlas_path = Path(atlas_path)

    if out_path is None:
        out_path = Path("rfu_out.tsv")
    out_path = Path(out_path)

    in_path = out_path.with_suffix(".in.tsv")
    features.loc[:, ["cell_id", "cdr3aa", "trbv"]].to_csv(in_path, sep="\t", index=False)

    cmd = [
        "Rscript",
        str(rscript_path),
        "--in",
        str(in_path),
        "--atlas",
        str(atlas_path),
        "--out",
        str(out_path),
    ] + list(extra_args)

    proc = subprocess.run(cmd, capture_output=True, text=True)

    if proc.returncode != 0:
        raise RuntimeError(
            "RFU Rscript failed.\n"
            f"cmd: {' '.join(cmd)}\n"
            f"returncode: {proc.returncode}\n"
            f"stdout:\n{proc.stdout}\n"
            f"stderr:\n{proc.stderr}\n"
        )

    df = pd.read_csv(out_path, sep="\t")
    return RFURunResult(df=df, stdout=proc.stdout, stderr=proc.stderr, returncode=proc.returncode)