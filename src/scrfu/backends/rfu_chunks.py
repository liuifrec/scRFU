from __future__ import annotations

import hashlib
import json
import os
import time
from collections.abc import Callable, Mapping, Sequence
from concurrent.futures import Future, ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

from ..rfu import RFURunResult, file_sha256

CHUNK_SCHEMA_VERSION = 1
WRAPPER_SCHEMA_VERSION = "headerless-stable-id-v1"
QUERY_COLUMNS = ("unique_sequence_id", "cdr3aa", "trbv")
REQUIRED_OUTPUT_COLUMNS = (
    "unique_sequence_id",
    "rfu_id",
    "rfu_label",
    "rfu_score",
    "pass_thr",
    "upstream_n_miss",
)


class RFUChunkError(RuntimeError):
    """Raised when a restartable RFU chunk cannot be completed safely."""

    def __init__(
        self,
        message: str,
        *,
        returncode: int = 1,
        stdout: str = "",
        stderr: str = "",
    ) -> None:
        super().__init__(message)
        self.returncode = returncode
        self.stdout = stdout
        self.stderr = stderr


@dataclass(frozen=True)
class RFURunIdentity:
    mode: str
    threshold: float
    deduplication_key: str | None
    eligibility_rule: str
    rfu_r_sha256: str
    trimer_rdata_sha256: str
    km5000_rdata_sha256: str
    wrapper_r_sha256: str
    scrfu_version: str
    chunk_size: int
    extra_r_args: tuple[str, ...] = ()
    wrapper_schema_version: str = WRAPPER_SCHEMA_VERSION

    def as_dict(self) -> dict[str, Any]:
        return {
            "mode": self.mode,
            "threshold": self.threshold,
            "deduplication_key": self.deduplication_key,
            "eligibility_rule": self.eligibility_rule,
            "rfu_r_sha256": self.rfu_r_sha256,
            "trimer_rdata_sha256": self.trimer_rdata_sha256,
            "km5000_rdata_sha256": self.km5000_rdata_sha256,
            "wrapper_r_sha256": self.wrapper_r_sha256,
            "scrfu_version": self.scrfu_version,
            "chunk_size": self.chunk_size,
            "extra_r_args": list(self.extra_r_args),
            "wrapper_schema_version": self.wrapper_schema_version,
        }


@dataclass(frozen=True)
class RFUChunkSpec:
    run_id: str
    chunk_id: str
    chunk_index: int
    start_offset: int
    end_offset_exclusive: int
    expected_row_count: int
    ordered_unique_sequence_ids: tuple[str, ...]
    ordered_identifier_sha256: str
    input_sha256: str


ChunkRunner = Callable[[pd.DataFrame, Path, Path, Path], RFURunResult]


def _run_chunk_task(
    runner: ChunkRunner,
    chunk: pd.DataFrame,
    input_path: Path,
    pending_output: Path,
    chunk_dir: Path,
) -> RFURunResult:
    """Top-level process-pool target; all paths are private to one chunk."""
    return runner(chunk, input_path, pending_output, chunk_dir)


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _canonical_json_bytes(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode(
        "utf-8"
    )


def _sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def _query_tsv_bytes(queries: pd.DataFrame) -> bytes:
    return queries.loc[:, list(QUERY_COLUMNS)].to_csv(sep="\t", index=False).encode("utf-8")


def _atomic_write_bytes(path: Path, value: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_bytes(value)
    os.replace(temporary, path)


def _atomic_write_text(path: Path, value: str) -> None:
    _atomic_write_bytes(path, value.encode("utf-8"))


def _atomic_write_json(path: Path, value: Mapping[str, Any]) -> None:
    _atomic_write_text(path, json.dumps(value, indent=2, sort_keys=True) + "\n")


def _optional_max_rss_kb() -> int | float | None:
    try:
        import resource

        return resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    except (ImportError, AttributeError, OSError, ValueError):
        return None


def deterministic_run_id(queries: pd.DataFrame, identity: RFURunIdentity) -> str:
    """Hash ordered scientific queries and execution-critical backend identity."""
    payload = {
        "ordered_unique_sequence_ids": queries["unique_sequence_id"].astype(str).tolist(),
        "ordered_cdr3aa": queries["cdr3aa"].astype(str).tolist(),
        **identity.as_dict(),
    }
    return _sha256_bytes(_canonical_json_bytes(payload))


def build_chunk_specs(
    queries: pd.DataFrame,
    *,
    run_id: str,
    chunk_size: int,
) -> list[RFUChunkSpec]:
    specs: list[RFUChunkSpec] = []
    for chunk_index, start in enumerate(range(0, len(queries), chunk_size)):
        end = min(start + chunk_size, len(queries))
        chunk = queries.iloc[start:end]
        identifiers = tuple(chunk["unique_sequence_id"].astype(str))
        identifier_sha = _sha256_bytes(_canonical_json_bytes(list(identifiers)))
        input_sha = _sha256_bytes(_query_tsv_bytes(chunk))
        chunk_id = f"{run_id[:12]}-{chunk_index:05d}-{input_sha[:12]}"
        specs.append(
            RFUChunkSpec(
                run_id=run_id,
                chunk_id=chunk_id,
                chunk_index=chunk_index,
                start_offset=start,
                end_offset_exclusive=end,
                expected_row_count=end - start,
                ordered_unique_sequence_ids=identifiers,
                ordered_identifier_sha256=identifier_sha,
                input_sha256=input_sha,
            )
        )
    return specs


def _chunk_manifest_base(
    spec: RFUChunkSpec,
    identity: RFURunIdentity,
) -> dict[str, Any]:
    return {
        "schema_version": CHUNK_SCHEMA_VERSION,
        "run_id": spec.run_id,
        "chunk_id": spec.chunk_id,
        "chunk_index": spec.chunk_index,
        "start_offset": spec.start_offset,
        "end_offset_exclusive": spec.end_offset_exclusive,
        "expected_row_count": spec.expected_row_count,
        "ordered_unique_sequence_ids": list(spec.ordered_unique_sequence_ids),
        "ordered_identifier_sha256": spec.ordered_identifier_sha256,
        "input_sha256": spec.input_sha256,
        "wrapper_sha256": identity.wrapper_r_sha256,
        "rfu_r_sha256": identity.rfu_r_sha256,
        "trimer_rdata_sha256": identity.trimer_rdata_sha256,
        "km5000_rdata_sha256": identity.km5000_rdata_sha256,
        "backend_mode": identity.mode,
        "threshold": identity.threshold,
        "eligibility_rule": identity.eligibility_rule,
        "deduplication_key": identity.deduplication_key,
        "scrfu_version": identity.scrfu_version,
        "wrapper_schema_version": identity.wrapper_schema_version,
    }


def _validate_output_table(
    output: pd.DataFrame,
    expected_identifiers: Sequence[str],
) -> str | None:
    missing = set(REQUIRED_OUTPUT_COLUMNS).difference(output.columns)
    if missing:
        return f"output missing required columns: {sorted(missing)}"
    if len(output) != len(expected_identifiers):
        return f"output row count mismatch: {len(output)} != {len(expected_identifiers)}"
    if output["unique_sequence_id"].duplicated().any():
        return "output contains duplicated unique_sequence_id values"
    actual = output["unique_sequence_id"].astype(str).tolist()
    if actual != list(expected_identifiers):
        missing_ids = sorted(set(expected_identifiers).difference(actual))
        unexpected_ids = sorted(set(actual).difference(expected_identifiers))
        return (
            f"output identifiers/order mismatch; missing={missing_ids}, unexpected={unexpected_ids}"
        )
    return None


def validate_cached_chunk(
    *,
    chunk_dir: Path,
    spec: RFUChunkSpec,
    identity: RFURunIdentity,
) -> tuple[pd.DataFrame | None, str]:
    manifest_path = chunk_dir / "manifest.json"
    output_path = chunk_dir / "output.tsv"
    input_path = chunk_dir / "input.tsv"
    if not manifest_path.exists():
        if output_path.exists() or input_path.exists():
            return None, "manifest missing while chunk artifacts exist"
        return None, "cache entry does not exist"
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        return None, f"manifest unreadable: {exc}"

    expected = _chunk_manifest_base(spec, identity)
    for key, value in expected.items():
        if manifest.get(key) != value:
            return None, f"manifest {key} mismatch"
    if manifest.get("status") != "complete":
        return None, f"manifest status is {manifest.get('status')!r}, not 'complete'"
    if not input_path.is_file():
        return None, "input.tsv missing"
    if file_sha256(input_path) != spec.input_sha256:
        return None, "input.tsv SHA256 mismatch"
    if not output_path.is_file():
        return None, "output.tsv missing"
    output_sha = file_sha256(output_path)
    if manifest.get("output_sha256") != output_sha:
        return None, "output.tsv SHA256 mismatch"
    try:
        output = pd.read_csv(output_path, sep="\t")
    except Exception as exc:
        return None, f"output.tsv unreadable: {exc}"
    output_error = _validate_output_table(output, spec.ordered_unique_sequence_ids)
    if output_error is not None:
        return None, output_error
    return output, "cache hit"


def _upstream_n_miss(output: pd.DataFrame) -> int:
    values = output["upstream_n_miss"].dropna().unique()
    if len(values) != 1:
        raise RFUChunkError("Chunk output has inconsistent upstream_n_miss values.")
    return int(values[0])


def execute_restartable_chunks(
    queries: pd.DataFrame,
    *,
    workdir: Path,
    identity: RFURunIdentity,
    runner: ChunkRunner,
    resume: bool,
    force_recompute: bool,
    max_workers: int = 1,
    executor: str = "process",
    run_context: Mapping[str, Any] | None = None,
) -> RFURunResult:
    if isinstance(max_workers, bool) or not isinstance(max_workers, int) or max_workers < 1:
        raise ValueError("max_workers must be a positive integer.")
    executor_name = str(executor).strip().lower()
    if executor_name not in {"process", "thread"}:
        raise ValueError("executor must be 'process' or 'thread'.")
    started = time.perf_counter()
    run_id = deterministic_run_id(queries, identity)
    specs = build_chunk_specs(queries, run_id=run_id, chunk_size=identity.chunk_size)
    run_dir = workdir.resolve() / "runs" / run_id
    chunks_dir = run_dir / "chunks"
    chunks_dir.mkdir(parents=True, exist_ok=True)
    run_manifest_path = run_dir / "run_manifest.json"

    counters = {
        "cache_hit_count": 0,
        "cache_miss_count": 0,
        "cache_invalidated_count": 0,
        "forced_recompute_count": 0,
        "recomputed_chunk_count": 0,
        "failed_chunk_count": 0,
    }
    cache_events: list[dict[str, Any]] = []
    outputs: dict[int, pd.DataFrame] = {}
    stdout_parts: dict[int, str] = {}
    stderr_parts: dict[int, str] = {}

    run_manifest: dict[str, Any] = {
        "schema_version": CHUNK_SCHEMA_VERSION,
        "wrapper_schema_version": identity.wrapper_schema_version,
        "run_id": run_id,
        "status": "running",
        "created_at": _utc_now(),
        "query_count": len(queries),
        "unique_query_count": len(queries),
        "chunk_size": identity.chunk_size,
        "chunk_count": len(specs),
        "max_workers": max_workers,
        "executor": executor_name,
        "identity": identity.as_dict(),
        "chunks": [spec.chunk_id for spec in specs],
        "run_manifest_path": str(run_manifest_path),
        **dict(run_context or {}),
    }
    _atomic_write_json(run_manifest_path, run_manifest)

    pending: list[dict[str, Any]] = []
    for spec in specs:
        chunk = queries.iloc[spec.start_offset : spec.end_offset_exclusive].reset_index(drop=True)
        chunk_dir = chunks_dir / f"chunk_{spec.chunk_index:05d}"
        chunk_dir.mkdir(parents=True, exist_ok=True)
        manifest_path = chunk_dir / "manifest.json"
        output_path = chunk_dir / "output.tsv"
        input_path = chunk_dir / "input.tsv"
        pending_output = chunk_dir / "output.tsv.pending"
        cache_artifacts_exist = (
            manifest_path.exists() or output_path.exists() or input_path.exists()
        )

        cached: pd.DataFrame | None = None
        reason = "cache entry does not exist"
        event = "cache_miss"
        if force_recompute:
            event = "forced_recomputation"
            reason = "force_recompute=True"
            counters["forced_recompute_count"] += 1
        elif resume:
            cached, reason = validate_cached_chunk(
                chunk_dir=chunk_dir,
                spec=spec,
                identity=identity,
            )
            if cached is not None:
                event = "cache_hit"
                counters["cache_hit_count"] += 1
            elif cache_artifacts_exist:
                event = "cache_invalidated"
                counters["cache_invalidated_count"] += 1
            else:
                counters["cache_miss_count"] += 1
        else:
            event = "resume_disabled_recomputation"
            reason = "resume=False"
            if not cache_artifacts_exist:
                counters["cache_miss_count"] += 1

        cache_events.append(
            {
                "chunk_id": spec.chunk_id,
                "chunk_index": spec.chunk_index,
                "event": event,
                "reason": reason,
            }
        )
        if cached is not None:
            outputs[spec.chunk_index] = cached
            continue

        counters["recomputed_chunk_count"] += 1
        if manifest_path.exists():
            try:
                previous = manifest_path.read_text(encoding="utf-8")
                _atomic_write_text(chunk_dir / "manifest.previous.json", previous)
            except OSError:
                pass
        _atomic_write_bytes(input_path, _query_tsv_bytes(chunk))
        manifest = {
            **_chunk_manifest_base(spec, identity),
            "status": "pending",
            "cache_event": event,
            "cache_invalidation_reason": reason if event == "cache_invalidated" else None,
            "start_timestamp": None,
            "completion_timestamp": None,
            "elapsed_seconds": None,
            "exit_code": None,
            "output_sha256": None,
            "process_max_rss_kb": None,
        }
        _atomic_write_json(manifest_path, manifest)
        chunk_started = time.perf_counter()
        manifest["status"] = "running"
        manifest["start_timestamp"] = _utc_now()
        _atomic_write_json(manifest_path, manifest)
        pending.append(
            {
                "spec": spec,
                "chunk": chunk,
                "chunk_dir": chunk_dir,
                "input_path": input_path,
                "pending_output": pending_output,
                "output_path": output_path,
                "manifest_path": manifest_path,
                "manifest": manifest,
                "started": chunk_started,
            }
        )

    def complete(task: dict[str, Any], result: RFURunResult) -> None:
        spec: RFUChunkSpec = task["spec"]
        if result.returncode != 0:
            raise RFUChunkError(
                "RFU runner returned a non-zero exit code.",
                returncode=result.returncode,
                stdout=result.stdout,
                stderr=result.stderr,
            )
        output_error = _validate_output_table(result.df, spec.ordered_unique_sequence_ids)
        if output_error is not None:
            raise RFUChunkError(output_error, stdout=result.stdout, stderr=result.stderr)
        _atomic_write_text(task["chunk_dir"] / "stdout.log", result.stdout or "")
        _atomic_write_text(task["chunk_dir"] / "stderr.log", result.stderr or "")
        _atomic_write_bytes(
            task["pending_output"], result.df.to_csv(sep="\t", index=False).encode()
        )
        os.replace(task["pending_output"], task["output_path"])
        manifest = task["manifest"]
        manifest.update(
            {
                "status": "complete",
                "completion_timestamp": _utc_now(),
                "elapsed_seconds": time.perf_counter() - task["started"],
                "exit_code": 0,
                "output_sha256": file_sha256(task["output_path"]),
                "process_max_rss_kb": _optional_max_rss_kb(),
            }
        )
        _atomic_write_json(task["manifest_path"], manifest)
        outputs[spec.chunk_index] = result.df
        stdout_parts[spec.chunk_index] = result.stdout
        stderr_parts[spec.chunk_index] = result.stderr

    def fail(task: dict[str, Any], exc: Exception) -> RFUChunkError:
        spec: RFUChunkSpec = task["spec"]
        counters["failed_chunk_count"] += 1
        returncode = int(getattr(exc, "returncode", 1))
        stdout = str(getattr(exc, "stdout", ""))
        stderr = str(getattr(exc, "stderr", ""))
        _atomic_write_text(task["chunk_dir"] / "stdout.log", stdout)
        _atomic_write_text(task["chunk_dir"] / "stderr.log", stderr)
        manifest = task["manifest"]
        manifest.update(
            {
                "status": "failed",
                "completion_timestamp": _utc_now(),
                "elapsed_seconds": time.perf_counter() - task["started"],
                "exit_code": returncode,
                "failure": str(exc),
                "process_max_rss_kb": _optional_max_rss_kb(),
            }
        )
        _atomic_write_json(task["manifest_path"], manifest)
        run_manifest.update(
            {
                "status": "failed",
                "completion_timestamp": _utc_now(),
                "total_elapsed_seconds": time.perf_counter() - started,
                "failed_chunk_id": spec.chunk_id,
                "completed_chunk_count": len(outputs),
                "reused_chunk_count": counters["cache_hit_count"],
                "invalidated_chunk_count": counters["cache_invalidated_count"],
                "cache_events": cache_events,
                **counters,
            }
        )
        _atomic_write_json(run_manifest_path, run_manifest)
        return RFUChunkError(
            f"RFU chunk {spec.chunk_id} (index {spec.chunk_index}) failed: {exc}",
            returncode=returncode,
            stdout=stdout,
            stderr=stderr,
        )

    if max_workers == 1:
        for task in pending:
            try:
                result = _run_chunk_task(
                    runner,
                    task["chunk"],
                    task["input_path"],
                    task["pending_output"],
                    task["chunk_dir"],
                )
                complete(task, result)
            except Exception as exc:
                raise fail(task, exc) from exc
    elif pending:
        pool_class = ProcessPoolExecutor if executor_name == "process" else ThreadPoolExecutor
        pool = pool_class(max_workers=max_workers)
        futures: dict[Future[RFURunResult], dict[str, Any]] = {
            pool.submit(
                _run_chunk_task,
                runner,
                task["chunk"],
                task["input_path"],
                task["pending_output"],
                task["chunk_dir"],
            ): task
            for task in pending
        }
        failure: tuple[dict[str, Any], Exception] | None = None
        try:
            for future in as_completed(futures):
                task = futures[future]
                try:
                    complete(task, future.result())
                except Exception as exc:
                    failure = (task, exc)
                    for other in futures:
                        if other is not future:
                            other.cancel()
                    break
        finally:
            pool.shutdown(wait=True, cancel_futures=True)
        if failure is not None:
            task, exc = failure
            raise fail(task, exc) from exc

    ordered_outputs = [outputs[index] for index in sorted(outputs)]
    assignments = (
        pd.concat(ordered_outputs, ignore_index=True) if ordered_outputs else pd.DataFrame()
    )
    global_error = _validate_output_table(
        assignments,
        queries["unique_sequence_id"].astype(str).tolist(),
    )
    if global_error is not None:
        raise RFUChunkError(f"Concatenated RFU output validation failed: {global_error}")
    total_elapsed = time.perf_counter() - started
    upstream_n_total = sum(_upstream_n_miss(output) for output in ordered_outputs)
    run_manifest.update(
        {
            "status": "complete",
            "completion_timestamp": _utc_now(),
            "total_elapsed_seconds": total_elapsed,
            "completed_chunk_count": len(specs),
            "reused_chunk_count": counters["cache_hit_count"],
            "invalidated_chunk_count": counters["cache_invalidated_count"],
            "cache_events": cache_events,
            **counters,
            "process_max_rss_kb": _optional_max_rss_kb(),
        }
    )
    _atomic_write_json(run_manifest_path, run_manifest)
    metadata = {
        "run_id": run_id,
        "chunking_enabled": True,
        "chunk_size": identity.chunk_size,
        "chunk_count": len(specs),
        "max_workers": max_workers,
        "executor": executor_name,
        "completed_chunk_count": len(specs),
        "reused_chunk_count": counters["cache_hit_count"],
        "recomputed_chunk_count": counters["recomputed_chunk_count"],
        "invalidated_chunk_count": counters["cache_invalidated_count"],
        "failed_chunk_count": counters["failed_chunk_count"],
        "cache_hit_count": counters["cache_hit_count"],
        "cache_miss_count": counters["cache_miss_count"],
        "forced_recompute_count": counters["forced_recompute_count"],
        "total_elapsed_seconds": total_elapsed,
        "run_manifest_path": str(run_manifest_path),
        "process_max_rss_kb": _optional_max_rss_kb(),
        "upstream_threshold_miss_count": upstream_n_total,
    }
    return RFURunResult(
        df=assignments,
        stdout="".join(stdout_parts[index] for index in sorted(stdout_parts)),
        stderr="".join(stderr_parts[index] for index in sorted(stderr_parts)),
        returncode=0,
        metadata=metadata,
    )


__all__ = [
    "CHUNK_SCHEMA_VERSION",
    "REQUIRED_OUTPUT_COLUMNS",
    "RFUChunkError",
    "RFUChunkSpec",
    "RFURunIdentity",
    "WRAPPER_SCHEMA_VERSION",
    "build_chunk_specs",
    "deterministic_run_id",
    "execute_restartable_chunks",
    "validate_cached_chunk",
]
