from __future__ import annotations

import os
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence, Union

import pandas as pd

from ..rfu import RFURunResult, file_sha256

PathLike = Union[str, Path]


@dataclass(frozen=True)
class RFURepoPaths:
    """
    Represents a user-provided checkout of the upstream RFU repository
    (https://github.com/s175573/RFU).

    We do NOT vendor RFU code or Rdata into scrfu. Users must supply rfu_dir.

    Required files inside rfu_dir:
      - RFU.R
      - trimerMDSfit_small.Rdata
      - km5000noMax.Rdata
    """

    rfu_dir: Path
    rfu_r: Path
    trimer_rdata: Path
    km5000_rdata: Path

    @staticmethod
    def from_dir(rfu_dir: PathLike) -> "RFURepoPaths":
        d = Path(rfu_dir).expanduser().resolve()
        rfu_r = d / "RFU.R"
        trimer = d / "trimerMDSfit_small.Rdata"
        km5000 = d / "km5000noMax.Rdata"

        missing = [p.name for p in [rfu_r, trimer, km5000] if not p.exists()]
        if missing:
            raise FileNotFoundError(
                "RFU_DIR is missing required files: "
                + ", ".join(missing)
                + f"\nRFU_DIR={d}\n"
                "Expected: RFU.R, trimerMDSfit_small.Rdata, km5000noMax.Rdata"
            )

        return RFURepoPaths(rfu_dir=d, rfu_r=rfu_r, trimer_rdata=trimer, km5000_rdata=km5000)


class RFURepoBackend:
    """
    Backend that calls the *upstream RFU repo code* via a thin wrapper R script.

    This backend expects:
      - features: DataFrame with columns [cell_id, cdr3aa, trbv]
    Produces:
      - output: DataFrame with columns [cell_id, rfu_label, rfu_score] (and possibly more)

    The wrapper R script is scrfu/r/run_rfu_repo.R (in this repository).
    """

    def __init__(
        self,
        *,
        rfu_dir: PathLike,
        wrapper_r_path: Optional[PathLike] = None,
        rscript_bin: str = "Rscript",
    ) -> None:
        self.paths = RFURepoPaths.from_dir(rfu_dir)

        if wrapper_r_path is None:
            # default: repo-relative path if user calls from installed package
            # users typically pass explicit path "r/run_rfu_repo.R" from their checkout
            wrapper_r_path = Path("r") / "run_rfu_repo.R"
        self.wrapper_r_path = Path(wrapper_r_path).expanduser().resolve()
        if not self.wrapper_r_path.exists():
            raise FileNotFoundError(
                f"Wrapper R script not found: {self.wrapper_r_path}\n"
                "Expected r/run_rfu_repo.R in your scrfu checkout."
            )

        self.rscript_bin = rscript_bin

    def run(
        self,
        features: pd.DataFrame,
        *,
        extra_args: Optional[Sequence[str]] = None,
        workdir: Optional[PathLike] = None,
    ) -> RFURunResult:
        """
        Run RFU assignment through the upstream repo code.

        Parameters
        ----------
        features
            DataFrame with columns: cell_id, cdr3aa, trbv
        extra_args
            Extra args forwarded to the wrapper R script (optional).
        workdir
            Optional directory to write temp files. If None, uses a temp dir.
        """
        required = {"cell_id", "cdr3aa", "trbv"}
        if not required.issubset(features.columns):
            raise ValueError(f"features must contain {required}, got {list(features.columns)}")

        extra_args = list(extra_args or [])

        # Work directory
        if workdir is None:
            td = tempfile.TemporaryDirectory(prefix="scrfu_rfu_")
            work_path = Path(td.name)
        else:
            work_path = Path(workdir).expanduser().resolve()
            work_path.mkdir(parents=True, exist_ok=True)
            td = None  # no tempdir cleanup

        in_tsv = work_path / "rfu_in.tsv"
        out_tsv = work_path / "rfu_out.tsv"

        features.loc[:, ["cell_id", "cdr3aa", "trbv"]].to_csv(in_tsv, sep="\t", index=False)

        cmd = [
            self.rscript_bin,
            str(self.wrapper_r_path),
            "--in",
            str(in_tsv),
            "--out",
            str(out_tsv),
            "--rfu_dir",
            str(self.paths.rfu_dir),
        ] + extra_args

        # Use the existing scrfu.run_rfu_rscript pattern but specialized here:
        import subprocess

        proc = subprocess.run(cmd, capture_output=True, text=True)

        # Cleanup temp dir if we created it
        try:
            if proc.returncode != 0:
                raise RuntimeError(
                    "RFURepoBackend failed.\n"
                    f"cmd: {' '.join(cmd)}\n"
                    f"returncode: {proc.returncode}\n"
                    f"stdout:\n{proc.stdout}\n"
                    f"stderr:\n{proc.stderr}\n"
                )

            if not out_tsv.exists():
                raise RuntimeError(
                    f"RFU wrapper did not produce output file: {out_tsv}\n"
                    f"stdout:\n{proc.stdout}\n"
                    f"stderr:\n{proc.stderr}\n"
                )

            df = pd.read_csv(out_tsv, sep="\t")
            return RFURunResult(df=df, stdout=proc.stdout, stderr=proc.stderr, returncode=proc.returncode)
        finally:
            if td is not None:
                td.cleanup()

    def provenance_dict(self) -> dict:
        """
        Return provenance info suitable for storing in adata.uns["scrfu"].
        """
        d = self.paths
        prov = {
            "backend": "rfu_repo",
            "rfu_dir": str(d.rfu_dir),
            "rfu_r_sha256": file_sha256(d.rfu_r),
            "trimer_rdata_sha256": file_sha256(d.trimer_rdata),
            "km5000_rdata_sha256": file_sha256(d.km5000_rdata),
            "wrapper_r_path": str(self.wrapper_r_path),
            "wrapper_r_sha256": file_sha256(self.wrapper_r_path) if self.wrapper_r_path.exists() else None,
            "rscript_bin": self.rscript_bin,
        }
        # Also store env hint (optional)
        prov["RFU_DIR_env"] = os.environ.get("RFU_DIR")
        return prov