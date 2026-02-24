from __future__ import annotations

import os
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence, Union

import pandas as pd

from ..rfu import RFURunResult, file_sha256

PathLike = Union[str, Path]


# ---------------------------------------------------------------------
# RFU upstream repository path validation
# ---------------------------------------------------------------------
@dataclass(frozen=True)
class RFURepoPaths:
    """
    Represents a user-provided checkout of the upstream RFU repository:
        https://github.com/s175573/RFU

    We do NOT vendor RFU code or Rdata into scrfu.
    Users must supply rfu_dir.
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

        missing = [p.name for p in (rfu_r, trimer, km5000) if not p.exists()]
        if missing:
            raise FileNotFoundError(
                "RFU_DIR is missing required files: "
                + ", ".join(missing)
                + f"\nRFU_DIR={d}\n"
                "Expected files:\n"
                "  - RFU.R\n"
                "  - trimerMDSfit_small.Rdata\n"
                "  - km5000noMax.Rdata"
            )

        return RFURepoPaths(
            rfu_dir=d,
            rfu_r=rfu_r,
            trimer_rdata=trimer,
            km5000_rdata=km5000,
        )


# ---------------------------------------------------------------------
# RFU Backend
# ---------------------------------------------------------------------
class RFURepoBackend:
    """
    Backend that calls the upstream RFU repo via scrfu's wrapper R script.

    Input:
        features DataFrame with columns:
            - cell_id
            - cdr3aa
            - trbv

    Output:
        DataFrame returned by R wrapper (must contain at least):
            - cell_id
            - rfu_label
            - rfu_score
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
            wrapper_r_path = Path("r") / "run_rfu_repo.R"

        self.wrapper_r_path = Path(wrapper_r_path).expanduser().resolve()

        if not self.wrapper_r_path.exists():
            raise FileNotFoundError(
                f"Wrapper R script not found:\n  {self.wrapper_r_path}\n"
                "Expected r/run_rfu_repo.R inside scrfu repo."
            )

        self.rscript_bin = rscript_bin

    # -----------------------------------------------------------------
    # Main execution
    # -----------------------------------------------------------------
    def run(
        self,
        features: pd.DataFrame,
        *,
        extra_args: Optional[Sequence[str]] = None,
        workdir: Optional[PathLike] = None,
    ) -> RFURunResult:

        required = {"cell_id", "cdr3aa", "trbv"}
        if not required.issubset(features.columns):
            raise ValueError(
                f"features must contain columns {required}, "
                f"but got {list(features.columns)}"
            )

        extra_args = list(extra_args or [])

        # -------------------------------------------------------------
        # Create deterministic repo-local workdir
        # -------------------------------------------------------------
        if workdir is None:
            base = Path(".scrfu") / "rfu_repo_runs"
            base.mkdir(parents=True, exist_ok=True)

            work_path = Path(
                tempfile.mkdtemp(prefix="run_", dir=str(base))
            ).resolve()
        else:
            work_path = Path(workdir).expanduser().resolve()
            work_path.mkdir(parents=True, exist_ok=True)

        in_tsv = work_path / "rfu_in.tsv"
        out_tsv = work_path / "rfu_out.tsv"
        stdout_log = work_path / "stdout.log"
        stderr_log = work_path / "stderr.log"

        # Write features
        features.loc[:, ["cell_id", "cdr3aa", "trbv"]].to_csv(
            in_tsv, sep="\t", index=False
        )

        cmd = [
            self.rscript_bin,
            str(self.wrapper_r_path),
            "--in",
            str(in_tsv),
            "--out",
            str(out_tsv),
            "--rfu_dir",
            str(self.paths.rfu_dir),
            "--workdir",
            str(work_path),
        ] + extra_args

        # -------------------------------------------------------------
        # Execute R
        # -------------------------------------------------------------
        proc = subprocess.run(cmd, capture_output=True, text=True)

        stdout_log.write_text(proc.stdout or "")
        stderr_log.write_text(proc.stderr or "")

        if proc.returncode != 0:
            raise RuntimeError(
                "RFURepoBackend execution failed.\n"
                f"Workdir: {work_path}\n"
                f"Command: {' '.join(cmd)}\n"
                f"Return code: {proc.returncode}\n"
                f"stdout log: {stdout_log}\n"
                f"stderr log: {stderr_log}\n"
            )

        if not out_tsv.exists():
            raise RuntimeError(
                "RFU wrapper did not produce expected output file.\n"
                f"Workdir: {work_path}\n"
                f"stdout log: {stdout_log}\n"
                f"stderr log: {stderr_log}\n"
            )

        df = pd.read_csv(out_tsv, sep="\t")

        return RFURunResult(
            df=df,
            stdout=proc.stdout,
            stderr=proc.stderr,
            returncode=proc.returncode,
        )

    # -----------------------------------------------------------------
    # Provenance
    # -----------------------------------------------------------------
    def provenance_dict(self) -> dict:
        d = self.paths

        return {
            "backend": "rfu_repo",
            "rfu_dir": str(d.rfu_dir),
            "rfu_r_sha256": file_sha256(d.rfu_r),
            "trimer_rdata_sha256": file_sha256(d.trimer_rdata),
            "km5000_rdata_sha256": file_sha256(d.km5000_rdata),
            "wrapper_r_path": str(self.wrapper_r_path),
            "wrapper_r_sha256": file_sha256(self.wrapper_r_path),
            "rscript_bin": self.rscript_bin,
            "RFU_DIR_env": os.environ.get("RFU_DIR"),
        }