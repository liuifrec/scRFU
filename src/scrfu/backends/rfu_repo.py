from __future__ import annotations

import os
import re
import subprocess
import tempfile
import time
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Literal, cast

import pandas as pd

from .._version import __version__
from ..rfu import RFURunResult, file_sha256
from .rfu_chunks import (
    RFUChunkError,
    RFURunIdentity,
    deterministic_run_id,
    execute_restartable_chunks,
)

PathLike = str | Path
RFUBackendMode = Literal["standard", "map_aware"]

_REQUIRED_FILES = ("RFU.R", "trimerMDSfit_small.Rdata", "km5000noMax.Rdata")
_FUNCTION_PATTERN = re.compile(
    r"^[ \t]*(AssignRFUs|AssignRFUs_with_map|RFUbatch_with_maps)"
    r"[ \t]*(?:<-|=)[ \t]*function[ \t]*\(",
    re.MULTILINE,
)


class RFUConfigurationError(ValueError):
    """Raised when an RFU checkout cannot be resolved from configuration."""


class RFUCapabilityError(RuntimeError):
    """Raised when an RFU checkout lacks a capability required by a mode."""


@dataclass(frozen=True)
class RFUCapabilities:
    """Functions detected statically in a supplied ``RFU.R`` file."""

    assign_rfus: bool
    assign_rfus_with_map: bool
    rfu_batch_with_maps: bool

    @classmethod
    def from_rfu_r(cls, path: PathLike) -> RFUCapabilities:
        text = Path(path).read_text(encoding="utf-8", errors="replace")
        functions = set(_FUNCTION_PATTERN.findall(text))
        return cls(
            assign_rfus="AssignRFUs" in functions,
            assign_rfus_with_map="AssignRFUs_with_map" in functions,
            rfu_batch_with_maps="RFUbatch_with_maps" in functions,
        )

    def as_dict(self) -> dict[str, bool]:
        return {
            "AssignRFUs": self.assign_rfus,
            "AssignRFUs_with_map": self.assign_rfus_with_map,
            "RFUbatch_with_maps": self.rfu_batch_with_maps,
        }

    def require(self, mode: RFUBackendMode) -> None:
        if mode == "standard" and not self.assign_rfus:
            raise RFUCapabilityError(
                "Standard RFU support was requested, but AssignRFUs() was not found in RFU.R. "
                "Supply an official-compatible RFU checkout containing AssignRFUs()."
            )
        if mode == "map_aware" and not self.assign_rfus_with_map:
            raise RFUCapabilityError(
                "Map-aware RFU support was requested, but AssignRFUs_with_map() was not found "
                'in RFU.R. Use mode="standard" or supply a compatible RFU checkout.'
            )

    def require_batch_with_maps(self) -> None:
        if not self.rfu_batch_with_maps:
            raise RFUCapabilityError(
                "RFU batch-map support was requested, but RFUbatch_with_maps() was not found "
                "in RFU.R. Use scRFU's standard per-repertoire workflow or supply a compatible "
                "optional checkout."
            )


@dataclass(frozen=True)
class RFURepoPaths:
    """Resolved and validated paths for a user-provided RFU checkout."""

    rfu_dir: Path
    rfu_r: Path
    trimer_rdata: Path
    km5000_rdata: Path
    resolution_source: Literal["explicit_argument", "environment_variable"]

    @classmethod
    def resolve(
        cls,
        rfu_dir: PathLike | None = None,
        *,
        environ: Mapping[str, str] | None = None,
    ) -> RFURepoPaths:
        environment = os.environ if environ is None else environ
        if rfu_dir is not None:
            selected = rfu_dir
            source = "explicit_argument"
        elif environment.get("RFU_DIR"):
            selected = environment["RFU_DIR"]
            source = "environment_variable"
        else:
            raise RFUConfigurationError(
                "RFU backend directory is not configured. Pass rfu_dir='/path/to/RFU' or set "
                "the RFU_DIR environment variable."
            )
        return cls.from_dir(selected, resolution_source=source)

    @classmethod
    def from_dir(
        cls,
        rfu_dir: PathLike,
        *,
        resolution_source: Literal[
            "explicit_argument", "environment_variable"
        ] = "explicit_argument",
    ) -> RFURepoPaths:
        directory = Path(rfu_dir).expanduser().resolve()
        rfu_r = directory / _REQUIRED_FILES[0]
        trimer = directory / _REQUIRED_FILES[1]
        km5000 = directory / _REQUIRED_FILES[2]
        required_paths = (rfu_r, trimer, km5000)

        missing = [path.name for path in required_paths if not path.is_file()]
        if missing:
            raise FileNotFoundError(
                "RFU_DIR is missing required files: "
                + ", ".join(missing)
                + f"\nRFU_DIR={directory}\nExpected files:\n"
                + "\n".join(f"  - {name}" for name in _REQUIRED_FILES)
            )

        return cls(
            rfu_dir=directory,
            rfu_r=rfu_r,
            trimer_rdata=trimer,
            km5000_rdata=km5000,
            resolution_source=resolution_source,
        )


def _validate_mode(mode: str) -> RFUBackendMode:
    normalized = str(mode).strip().lower()
    if normalized not in {"standard", "map_aware"}:
        raise ValueError(f"Unknown RFU backend mode: {mode!r}. Supported: 'standard', 'map_aware'.")
    return cast(RFUBackendMode, normalized)


def _prepare_queries(
    features: pd.DataFrame,
    *,
    deduplicate: bool,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    required = {"cell_id", "cdr3aa", "trbv"}
    if not required.issubset(features.columns):
        raise ValueError(
            f"features must contain columns {required}, but got {list(features.columns)}"
        )

    original = features.copy().reset_index(drop=True)
    original["input_row_id"] = range(len(original))
    cdr3 = original["cdr3aa"].astype("string")
    eligible = cdr3.notna() & cdr3.str.startswith("C", na=False)
    original["eligibility_status"] = "ineligible_cdr3_not_starting_c"
    original.loc[eligible, "eligibility_status"] = "eligible"
    original["unique_sequence_id"] = pd.Series(pd.NA, index=original.index, dtype="string")

    eligible_rows = original.loc[eligible, ["input_row_id", "cdr3aa", "trbv"]].copy()
    if deduplicate:
        queries = eligible_rows.drop_duplicates(subset=["cdr3aa"], keep="first").copy()
    else:
        queries = eligible_rows.copy()
    queries = queries.reset_index(drop=True)
    queries["unique_sequence_id"] = [f"sequence_{index:08d}" for index in range(len(queries))]

    if deduplicate:
        sequence_lookup = queries.set_index("cdr3aa")["unique_sequence_id"]
        original.loc[eligible, "unique_sequence_id"] = (
            original.loc[eligible, "cdr3aa"].map(sequence_lookup).astype("string")
        )
    else:
        original.loc[eligible, "unique_sequence_id"] = queries["unique_sequence_id"].to_numpy()

    return original, queries.loc[:, ["unique_sequence_id", "cdr3aa", "trbv"]]


def _reconstruct_results(
    original: pd.DataFrame,
    assignments: pd.DataFrame,
) -> pd.DataFrame:
    expected = {"unique_sequence_id", "rfu_id", "rfu_label", "rfu_score", "pass_thr"}
    if not expected.issubset(assignments.columns):
        raise RuntimeError(
            "RFU wrapper output is missing required columns: "
            f"{sorted(expected.difference(assignments.columns))}"
        )
    if assignments["unique_sequence_id"].duplicated().any():
        raise RuntimeError("RFU wrapper returned duplicate unique_sequence_id values.")

    assignment_columns = [
        "unique_sequence_id",
        "rfu_id",
        "rfu_label",
        "rfu_score",
        "pass_thr",
    ]
    result = original.merge(
        assignments.loc[:, assignment_columns],
        on="unique_sequence_id",
        how="left",
        sort=False,
        validate="many_to_one",
    )
    result = result.sort_values("input_row_id", kind="stable").reset_index(drop=True)

    eligible = result["eligibility_status"].eq("eligible")
    returned_ids = set(assignments["unique_sequence_id"].astype(str))
    expected_ids = set(result.loc[eligible, "unique_sequence_id"].dropna().astype(str))
    if returned_ids != expected_ids:
        raise RuntimeError(
            "RFU wrapper did not return exactly one result for every eligible query."
        )
    result["pass_thr"] = result["pass_thr"].astype("boolean")
    result["rfu_status"] = "ineligible_cdr3_not_starting_c"
    assigned = eligible & result["rfu_id"].notna()
    result.loc[eligible & result["rfu_id"].isna(), "rfu_status"] = "upstream_unassigned"
    result.loc[assigned & result["pass_thr"].fillna(False), "rfu_status"] = (
        "assigned_threshold_pass"
    )
    result.loc[assigned & ~result["pass_thr"].fillna(False), "rfu_status"] = (
        "assigned_below_threshold"
    )
    return result


class RFURepoBackend:
    """Run a validated external RFU checkout through scRFU's R wrapper."""

    def __init__(
        self,
        *,
        rfu_dir: PathLike | None = None,
        mode: str = "standard",
        wrapper_r_path: PathLike | None = None,
        rscript_bin: str = "Rscript",
        environ: Mapping[str, str] | None = None,
    ) -> None:
        self.paths = RFURepoPaths.resolve(rfu_dir, environ=environ)
        self.mode = _validate_mode(mode)
        self.capabilities = RFUCapabilities.from_rfu_r(self.paths.rfu_r)
        self.capabilities.require(self.mode)

        if wrapper_r_path is None:
            wrapper_r_path = Path("r") / "run_rfu_repo.R"
        self.wrapper_r_path = Path(wrapper_r_path).expanduser().resolve()
        if not self.wrapper_r_path.is_file():
            raise FileNotFoundError(
                f"Wrapper R script not found:\n  {self.wrapper_r_path}\n"
                "Expected r/run_rfu_repo.R inside scrfu repo."
            )
        self.rscript_bin = rscript_bin

    def _run_query_file(
        self,
        queries: pd.DataFrame,
        input_path: Path,
        output_path: Path,
        execution_dir: Path,
        *,
        threshold: float,
        extra_args: Sequence[str],
    ) -> RFURunResult:
        cmd = [
            self.rscript_bin,
            str(self.wrapper_r_path),
            "--in",
            str(input_path),
            "--out",
            str(output_path),
            "--rfu_dir",
            str(self.paths.rfu_dir),
            "--workdir",
            str(execution_dir),
            "--mode",
            self.mode,
            "--threshold",
            str(float(threshold)),
        ] + list(extra_args)
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            raise RFUChunkError(
                "RFURepoBackend execution failed.\n"
                f"Workdir: {execution_dir}\nCommand: {' '.join(cmd)}\n"
                f"Return code: {proc.returncode}\nstdout:\n{proc.stdout}\nstderr:\n{proc.stderr}",
                returncode=proc.returncode,
                stdout=proc.stdout,
                stderr=proc.stderr,
            )
        if not output_path.is_file():
            raise RFUChunkError(
                f"RFU wrapper did not produce expected output file: {output_path}",
                stdout=proc.stdout,
                stderr=proc.stderr,
            )
        assignments = pd.read_csv(output_path, sep="\t")
        required = {
            "unique_sequence_id",
            "rfu_id",
            "rfu_label",
            "rfu_score",
            "pass_thr",
            "upstream_n_miss",
        }
        missing = required.difference(assignments.columns)
        if missing:
            raise RFUChunkError(
                f"RFU wrapper output is missing required columns: {sorted(missing)}",
                stdout=proc.stdout,
                stderr=proc.stderr,
            )
        expected_ids = queries["unique_sequence_id"].astype(str).tolist()
        actual_ids = assignments["unique_sequence_id"].astype(str).tolist()
        if assignments["unique_sequence_id"].duplicated().any() or actual_ids != expected_ids:
            raise RFUChunkError(
                "RFU wrapper output identifiers are duplicated, missing, or out of order.",
                stdout=proc.stdout,
                stderr=proc.stderr,
            )
        return RFURunResult(
            df=assignments,
            stdout=proc.stdout,
            stderr=proc.stderr,
            returncode=proc.returncode,
        )

    def run(
        self,
        features: pd.DataFrame,
        *,
        threshold: float = 0.6,
        deduplicate: bool = True,
        chunk_size: int | None = None,
        resume: bool = True,
        force_recompute: bool = False,
        extra_args: Sequence[str] | None = None,
        workdir: PathLike | None = None,
    ) -> RFURunResult:
        started = time.perf_counter()
        if isinstance(extra_args, str):
            raise TypeError("extra_args must be a sequence of strings, not a single string.")
        if chunk_size is not None and (
            isinstance(chunk_size, bool) or not isinstance(chunk_size, int) or chunk_size <= 0
        ):
            raise ValueError("chunk_size must be a positive integer or None.")
        original, queries = _prepare_queries(features, deduplicate=deduplicate)
        extra_args = list(extra_args or [])

        artifact_hashes = {
            "rfu_r_sha256": file_sha256(self.paths.rfu_r),
            "trimer_rdata_sha256": file_sha256(self.paths.trimer_rdata),
            "km5000_rdata_sha256": file_sha256(self.paths.km5000_rdata),
            "wrapper_r_sha256": file_sha256(self.wrapper_r_path),
        }
        identity = RFURunIdentity(
            mode=self.mode,
            threshold=float(threshold),
            deduplication_key="cdr3aa" if deduplicate else None,
            eligibility_rule="^C",
            rfu_r_sha256=artifact_hashes["rfu_r_sha256"],
            trimer_rdata_sha256=artifact_hashes["trimer_rdata_sha256"],
            km5000_rdata_sha256=artifact_hashes["km5000_rdata_sha256"],
            wrapper_r_sha256=artifact_hashes["wrapper_r_sha256"],
            scrfu_version=__version__,
            chunk_size=chunk_size or 0,
            extra_r_args=tuple(extra_args),
        )
        run_id = deterministic_run_id(queries, identity)

        metadata = {
            "run_id": run_id,
            "eligibility_rule": "^C",
            "deduplication_key": "cdr3aa" if deduplicate else None,
            "original_row_count": len(original),
            "eligible_row_count": int(original["eligibility_status"].eq("eligible").sum()),
            "unique_query_count": len(queries),
            "deduplication_ratio": len(queries)
            / len(original.loc[original["eligibility_status"].eq("eligible")])
            if original["eligibility_status"].eq("eligible").any()
            else None,
            "rfu_threshold": float(threshold),
            "chunking_enabled": chunk_size is not None,
            "chunk_size": chunk_size,
        }

        if queries.empty:
            assignments = pd.DataFrame(
                columns=["unique_sequence_id", "rfu_id", "rfu_label", "rfu_score", "pass_thr"]
            )
            reconstructed = _reconstruct_results(original, assignments)
            metadata.update(
                {
                    "reconstructed_output_row_count": len(reconstructed),
                    "upstream_threshold_miss_count": 0,
                    "reconstructed_threshold_miss_count": 0,
                    "chunk_count": 0,
                    "completed_chunk_count": 0,
                    "reused_chunk_count": 0,
                    "recomputed_chunk_count": 0,
                    "invalidated_chunk_count": 0,
                    "failed_chunk_count": 0,
                    "cache_hit_count": 0,
                    "cache_miss_count": 0,
                    "forced_recompute_count": 0,
                    "total_elapsed_seconds": time.perf_counter() - started,
                    "run_manifest_path": None,
                    "process_max_rss_kb": None,
                }
            )
            return RFURunResult(
                df=reconstructed,
                stdout="",
                stderr="",
                returncode=0,
                metadata=metadata,
            )

        if chunk_size is not None:
            work_path = (
                Path(workdir).expanduser().resolve()
                if workdir is not None
                else (Path(".scrfu") / "rfu_repo_runs").resolve()
            )
            work_path.mkdir(parents=True, exist_ok=True)

            def runner(
                chunk: pd.DataFrame,
                input_path: Path,
                pending_output: Path,
                chunk_dir: Path,
            ) -> RFURunResult:
                return self._run_query_file(
                    chunk,
                    input_path,
                    pending_output,
                    chunk_dir,
                    threshold=threshold,
                    extra_args=extra_args,
                )

            query_run = execute_restartable_chunks(
                queries,
                workdir=work_path,
                identity=identity,
                runner=runner,
                resume=resume,
                force_recompute=force_recompute,
                run_context={
                    "original_row_count": len(original),
                    "eligible_row_count": metadata["eligible_row_count"],
                    "reconstructed_output_row_count": len(original),
                    "deduplication_ratio": metadata["deduplication_ratio"],
                    "rfu_threshold": float(threshold),
                },
            )
            assignments = query_run.df
            upstream_n = int(query_run.metadata["upstream_threshold_miss_count"])
            metadata.update(query_run.metadata)
            stdout = query_run.stdout
            stderr = query_run.stderr
        else:
            if workdir is None:
                base = Path(".scrfu") / "rfu_repo_runs"
                base.mkdir(parents=True, exist_ok=True)
                work_path = Path(tempfile.mkdtemp(prefix="run_", dir=str(base))).resolve()
            else:
                work_path = Path(workdir).expanduser().resolve()
                work_path.mkdir(parents=True, exist_ok=True)
            in_tsv = work_path / "rfu_in.tsv"
            out_tsv = work_path / "rfu_out.tsv"
            queries.to_csv(in_tsv, sep="\t", index=False)
            query_run = self._run_query_file(
                queries,
                in_tsv,
                out_tsv,
                work_path,
                threshold=threshold,
                extra_args=extra_args,
            )
            (work_path / "stdout.log").write_text(query_run.stdout or "")
            (work_path / "stderr.log").write_text(query_run.stderr or "")
            assignments = query_run.df
            values = assignments["upstream_n_miss"].dropna().unique()
            if len(values) != 1:
                raise RuntimeError("RFU wrapper returned inconsistent upstream_n_miss values.")
            upstream_n = int(values[0])
            stdout = query_run.stdout
            stderr = query_run.stderr
            metadata.update(
                {
                    "chunk_count": 1,
                    "completed_chunk_count": 1,
                    "reused_chunk_count": 0,
                    "recomputed_chunk_count": 1,
                    "invalidated_chunk_count": 0,
                    "failed_chunk_count": 0,
                    "cache_hit_count": 0,
                    "cache_miss_count": 0,
                    "forced_recompute_count": 0,
                    "run_manifest_path": None,
                    "process_max_rss_kb": None,
                }
            )

        reconstructed = _reconstruct_results(original, assignments)
        eligible_result = reconstructed[reconstructed["eligibility_status"].eq("eligible")]
        reconstructed_n = int((~eligible_result["pass_thr"].fillna(False)).sum())
        metadata.update(
            {
                "reconstructed_output_row_count": len(reconstructed),
                "upstream_threshold_miss_count": upstream_n,
                "reconstructed_threshold_miss_count": reconstructed_n,
                "total_elapsed_seconds": time.perf_counter() - started,
            }
        )
        return RFURunResult(
            df=reconstructed,
            stdout=stdout,
            stderr=stderr,
            returncode=0,
            metadata=metadata,
        )

    def provenance_dict(self) -> dict:
        paths = self.paths
        return {
            "scrfu_version": __version__,
            "backend": "rfu_repo",
            "backend_mode": self.mode,
            "backend_resolution_source": paths.resolution_source,
            "rfu_dir": str(paths.rfu_dir),
            "resolved_backend_path": str(paths.rfu_dir),
            "rfu_r_sha256": file_sha256(paths.rfu_r),
            "trimer_rdata_sha256": file_sha256(paths.trimer_rdata),
            "km5000_rdata_sha256": file_sha256(paths.km5000_rdata),
            "detected_capabilities": self.capabilities.as_dict(),
            "wrapper_r_path": str(self.wrapper_r_path),
            "wrapper_r_sha256": file_sha256(self.wrapper_r_path),
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "rscript_bin": self.rscript_bin,
        }


__all__ = [
    "RFUBackendMode",
    "RFUCapabilities",
    "RFUCapabilityError",
    "RFUConfigurationError",
    "RFURepoBackend",
    "RFURepoPaths",
]
