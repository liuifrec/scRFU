from __future__ import annotations

import hashlib
import json
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from scrfu.backends.rfu_chunks import (
    RFURunIdentity,
    build_chunk_specs,
    deterministic_run_id,
    execute_restartable_chunks,
)
from scrfu.backends.rfu_repo import RFURepoBackend
from scrfu.rfu import RFURunResult


@dataclass(frozen=True)
class _PickleableChunkRunner:
    reverse_completion: bool = False
    fail_index: int | None = None

    def __call__(
        self,
        queries: pd.DataFrame,
        input_path: Path,
        output_path: Path,
        execution_dir: Path,
    ) -> RFURunResult:
        del input_path, output_path
        index = int(execution_dir.name.rsplit("_", maxsplit=1)[-1])
        if self.reverse_completion:
            time.sleep(max(0, 3 - index) * 0.01)
        if self.fail_index == index:
            from scrfu.backends.rfu_chunks import RFUChunkError

            raise RFUChunkError("parallel injected failure", returncode=19)
        scores = [0.7 + index * 0.01] * len(queries)
        output = pd.DataFrame(
            {
                "unique_sequence_id": queries["unique_sequence_id"].astype(str),
                "rfu_id": [index + 1] * len(queries),
                "rfu_label": [f"RFU{index + 1}"] * len(queries),
                "rfu_score": scores,
                "pass_thr": [True] * len(queries),
                "upstream_n_miss": [0] * len(queries),
            }
        )
        return RFURunResult(output, f"chunk {index}\n", "", 0)


def _backend(tmp_path: Path) -> RFURepoBackend:
    rfu_dir = tmp_path / "RFU"
    rfu_dir.mkdir(parents=True)
    (rfu_dir / "RFU.R").write_text("AssignRFUs <- function(ff) {}\n", encoding="utf-8")
    (rfu_dir / "trimerMDSfit_small.Rdata").write_bytes(b"trimer")
    (rfu_dir / "km5000noMax.Rdata").write_bytes(b"centers")
    wrapper = tmp_path / "wrapper.R"
    wrapper.write_text("#!/usr/bin/env Rscript\n", encoding="utf-8")
    return RFURepoBackend(rfu_dir=rfu_dir, wrapper_r_path=wrapper)


def _features() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "cell_id": ["c0", "c1", "c2", "c3", "c4", "c5", "c6"],
            "cdr3aa": ["CASSA", "CASST", "CASSA", "BAD", "CASSQ", "CASSR", "CASSN"],
            "trbv": ["TRBV1", "TRBV2", "TRBV3", "TRBV4", "TRBV5", "TRBV6", "TRBV7"],
        }
    )


def _install_fake_runner(
    backend: RFURepoBackend,
    monkeypatch: pytest.MonkeyPatch,
    calls: list[str],
    *,
    fail_once_at: int | None = None,
) -> None:
    failed = False

    def fake_runner(
        queries: pd.DataFrame,
        input_path: Path,
        output_path: Path,
        execution_dir: Path,
        *,
        threshold: float,
        extra_args: list[str],
    ) -> RFURunResult:
        nonlocal failed
        del input_path, output_path, extra_args
        calls.append(execution_dir.name)
        chunk_index = (
            int(execution_dir.name.rsplit("_", maxsplit=1)[-1])
            if "chunk_" in execution_dir.name
            else -1
        )
        if fail_once_at == chunk_index and not failed:
            failed = True
            from scrfu.backends.rfu_chunks import RFUChunkError

            raise RFUChunkError(
                "injected failure",
                returncode=17,
                stdout="before failure",
                stderr="synthetic failure",
            )

        scores = []
        ids = []
        for cdr3 in queries["cdr3aa"].astype(str):
            digest = hashlib.sha256(cdr3.encode()).hexdigest()
            ids.append(int(digest[:6], 16) % 5000 + 1)
            scores.append(0.4 + (int(digest[6:10], 16) % 500) / 1000)
        output = pd.DataFrame(
            {
                "unique_sequence_id": queries["unique_sequence_id"].astype(str),
                "rfu_id": ids,
                "rfu_label": [f"RFU{value}" for value in ids],
                "rfu_score": scores,
                "pass_thr": [score >= threshold for score in scores],
                "upstream_n_miss": [sum(score < threshold for score in scores)] * len(queries),
            }
        )
        return RFURunResult(output, "fake stdout\n", "", 0)

    monkeypatch.setattr(backend, "_run_query_file", fake_runner)


def _identity(**updates: Any) -> RFURunIdentity:
    values: dict[str, Any] = {
        "mode": "standard",
        "threshold": 0.6,
        "deduplication_key": "cdr3aa",
        "eligibility_rule": "^C",
        "rfu_r_sha256": "a" * 64,
        "trimer_rdata_sha256": "b" * 64,
        "km5000_rdata_sha256": "c" * 64,
        "wrapper_r_sha256": "d" * 64,
        "scrfu_version": "0.test",
        "chunk_size": 2,
    }
    values.update(updates)
    return RFURunIdentity(**values)


def _queries() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "unique_sequence_id": ["sequence_00000000", "sequence_00000001", "sequence_00000002"],
            "cdr3aa": ["CASSA", "CASST", "CASSQ"],
            "trbv": ["TRBV1", "TRBV2", "TRBV3"],
        }
    )


@pytest.mark.parametrize("chunk_size", [0, -1, 1.5, True])
def test_invalid_chunk_size_values(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, chunk_size: Any
) -> None:
    backend = _backend(tmp_path)
    _install_fake_runner(backend, monkeypatch, [])

    with pytest.raises(ValueError, match="positive integer"):
        backend.run(_features(), chunk_size=chunk_size, workdir=tmp_path / "work")


def test_run_and_chunk_ids_are_deterministic_and_final_chunk_is_partial() -> None:
    queries = _queries()
    identity = _identity()

    run_a = deterministic_run_id(queries, identity)
    run_b = deterministic_run_id(queries.copy(), identity)
    specs_a = build_chunk_specs(queries, run_id=run_a, chunk_size=2)
    specs_b = build_chunk_specs(queries, run_id=run_b, chunk_size=2)

    assert run_a == run_b
    assert [spec.chunk_id for spec in specs_a] == [spec.chunk_id for spec in specs_b]
    assert [spec.expected_row_count for spec in specs_a] == [2, 1]
    assert specs_a[-1].start_offset == 2
    assert specs_a[-1].end_offset_exclusive == 3


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("threshold", 0.7),
        ("mode", "map_aware"),
        ("wrapper_r_sha256", "e" * 64),
        ("rfu_r_sha256", "f" * 64),
        ("trimer_rdata_sha256", "1" * 64),
        ("km5000_rdata_sha256", "2" * 64),
    ],
)
def test_run_id_changes_for_execution_critical_inputs(field: str, value: Any) -> None:
    queries = _queries()

    assert deterministic_run_id(queries, _identity()) != deterministic_run_id(
        queries, _identity(**{field: value})
    )


def test_run_id_changes_when_ordered_query_input_changes() -> None:
    queries = _queries()
    changed = queries.copy()
    changed.loc[1, "cdr3aa"] = "CASS_CHANGED"

    assert deterministic_run_id(queries, _identity()) != deterministic_run_id(changed, _identity())


@pytest.mark.parametrize("chunk_size", [1, 2, 4])
def test_chunked_and_unchunked_outputs_are_identical(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, chunk_size: int
) -> None:
    backend = _backend(tmp_path)
    _install_fake_runner(backend, monkeypatch, [])

    unchunked = backend.run(_features(), workdir=tmp_path / "single")
    chunked = backend.run(
        _features(), chunk_size=chunk_size, workdir=tmp_path / f"chunks-{chunk_size}"
    )

    pd.testing.assert_frame_equal(unchunked.df, chunked.df, check_exact=False, atol=1e-12, rtol=0)
    assert chunked.df["input_row_id"].tolist() == list(range(len(_features())))
    assert chunked.df["cell_id"].tolist() == _features()["cell_id"].tolist()


def test_valid_cache_is_reused(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    backend = _backend(tmp_path)
    calls: list[str] = []
    _install_fake_runner(backend, monkeypatch, calls)
    workdir = tmp_path / "work"

    first = backend.run(_features(), chunk_size=2, workdir=workdir)
    first_call_count = len(calls)
    second = backend.run(_features(), chunk_size=2, workdir=workdir)

    assert len(calls) == first_call_count
    assert second.metadata["reused_chunk_count"] == first.metadata["chunk_count"]
    assert second.metadata["recomputed_chunk_count"] == 0
    run_manifest = json.loads(Path(second.metadata["run_manifest_path"]).read_text())
    assert run_manifest["original_row_count"] == len(_features())
    assert run_manifest["eligible_row_count"] == 6
    assert run_manifest["unique_query_count"] == 5
    assert run_manifest["reconstructed_output_row_count"] == len(_features())
    assert run_manifest["reused_chunk_count"] == first.metadata["chunk_count"]
    pd.testing.assert_frame_equal(first.df, second.df)


@pytest.mark.parametrize(
    ("kwargs", "counter"),
    [
        ({"resume": False}, "recomputed_chunk_count"),
        ({"force_recompute": True}, "forced_recompute_count"),
    ],
)
def test_resume_controls_recomputation(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    kwargs: dict[str, bool],
    counter: str,
) -> None:
    backend = _backend(tmp_path)
    calls: list[str] = []
    _install_fake_runner(backend, monkeypatch, calls)
    workdir = tmp_path / "work"
    first = backend.run(_features(), chunk_size=2, workdir=workdir)
    first_call_count = len(calls)

    second = backend.run(_features(), chunk_size=2, workdir=workdir, **kwargs)

    assert len(calls) == first_call_count + first.metadata["chunk_count"]
    assert second.metadata[counter] == first.metadata["chunk_count"]
    assert second.metadata["reused_chunk_count"] == 0


def _only_chunk_dir(run: RFURunResult) -> Path:
    manifest = Path(run.metadata["run_manifest_path"])
    chunks = sorted((manifest.parent / "chunks").iterdir())
    assert len(chunks) == 1
    return chunks[0]


@pytest.mark.parametrize(
    "corruption",
    [
        "hash",
        "row_count",
        "missing_column",
        "missing_identifier",
        "duplicate_identifier",
    ],
)
def test_corrupt_output_invalidates_cache(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    corruption: str,
) -> None:
    backend = _backend(tmp_path)
    calls: list[str] = []
    _install_fake_runner(backend, monkeypatch, calls)
    workdir = tmp_path / "work"
    first = backend.run(_features(), chunk_size=20, workdir=workdir)
    chunk_dir = _only_chunk_dir(first)
    output_path = chunk_dir / "output.tsv"
    output = pd.read_csv(output_path, sep="\t")
    if corruption == "hash":
        output_path.write_text(output_path.read_text() + "\n", encoding="utf-8")
    elif corruption == "row_count":
        output.iloc[:-1].to_csv(output_path, sep="\t", index=False)
        manifest = json.loads((chunk_dir / "manifest.json").read_text())
        manifest["output_sha256"] = hashlib.sha256(output_path.read_bytes()).hexdigest()
        (chunk_dir / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    elif corruption == "missing_column":
        output.drop(columns="rfu_score").to_csv(output_path, sep="\t", index=False)
        manifest = json.loads((chunk_dir / "manifest.json").read_text())
        manifest["output_sha256"] = hashlib.sha256(output_path.read_bytes()).hexdigest()
        (chunk_dir / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    elif corruption == "missing_identifier":
        output.loc[0, "unique_sequence_id"] = "missing"
        output.to_csv(output_path, sep="\t", index=False)
        manifest = json.loads((chunk_dir / "manifest.json").read_text())
        manifest["output_sha256"] = hashlib.sha256(output_path.read_bytes()).hexdigest()
        (chunk_dir / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    else:
        output.loc[1, "unique_sequence_id"] = output.loc[0, "unique_sequence_id"]
        output.to_csv(output_path, sep="\t", index=False)
        manifest = json.loads((chunk_dir / "manifest.json").read_text())
        manifest["output_sha256"] = hashlib.sha256(output_path.read_bytes()).hexdigest()
        (chunk_dir / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")

    second = backend.run(_features(), chunk_size=20, workdir=workdir)

    assert second.metadata["invalidated_chunk_count"] == 1
    assert second.metadata["recomputed_chunk_count"] == 1


def test_changed_input_gets_a_new_run_and_does_not_reuse_old_chunks(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    backend = _backend(tmp_path)
    calls: list[str] = []
    _install_fake_runner(backend, monkeypatch, calls)
    workdir = tmp_path / "work"
    first = backend.run(_features(), chunk_size=2, workdir=workdir)
    first_call_count = len(calls)
    changed = _features()
    changed.loc[1, "cdr3aa"] = "CASS_CHANGED"

    second = backend.run(changed, chunk_size=2, workdir=workdir)

    assert second.metadata["run_id"] != first.metadata["run_id"]
    assert second.metadata["reused_chunk_count"] == 0
    assert second.metadata["cache_miss_count"] == second.metadata["chunk_count"]
    assert len(calls) == first_call_count + second.metadata["chunk_count"]


@pytest.mark.parametrize(
    "manifest_key",
    [
        "threshold",
        "backend_mode",
        "wrapper_sha256",
        "rfu_r_sha256",
        "trimer_rdata_sha256",
        "km5000_rdata_sha256",
    ],
)
def test_manifest_identity_mismatch_invalidates_cache(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    manifest_key: str,
) -> None:
    backend = _backend(tmp_path)
    _install_fake_runner(backend, monkeypatch, [])
    workdir = tmp_path / "work"
    first = backend.run(_features(), chunk_size=20, workdir=workdir)
    chunk_dir = _only_chunk_dir(first)
    manifest_path = chunk_dir / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest[manifest_key] = "changed"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    second = backend.run(_features(), chunk_size=20, workdir=workdir)

    assert second.metadata["invalidated_chunk_count"] == 1
    assert second.metadata["recomputed_chunk_count"] == 1


def test_failed_middle_chunk_resumes_without_repeating_completed_chunk(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    backend = _backend(tmp_path)
    calls: list[str] = []
    _install_fake_runner(backend, monkeypatch, calls, fail_once_at=1)
    workdir = tmp_path / "interrupted"

    with pytest.raises(RuntimeError, match=r"index 1.*injected failure"):
        backend.run(_features(), chunk_size=2, workdir=workdir)
    assert calls == ["chunk_00000", "chunk_00001"]

    resumed = backend.run(_features(), chunk_size=2, workdir=workdir, resume=True)
    assert calls.count("chunk_00000") == 1
    assert resumed.metadata["reused_chunk_count"] == 1
    assert resumed.metadata["completed_chunk_count"] == 3

    clean_backend = _backend(tmp_path / "clean-root")
    _install_fake_runner(clean_backend, monkeypatch, [])
    uninterrupted = clean_backend.run(_features(), chunk_size=2, workdir=tmp_path / "uninterrupted")
    pd.testing.assert_frame_equal(resumed.df, uninterrupted.df)


@pytest.mark.parametrize("executor", ["thread", "process"])
def test_parallel_chunks_match_serial_and_preserve_query_order(
    tmp_path: Path, executor: str
) -> None:
    queries = pd.concat(
        [_queries(), _queries().assign(unique_sequence_id=lambda x: x["unique_sequence_id"] + "x")],
        ignore_index=True,
    )
    identity = _identity(chunk_size=2)
    serial = execute_restartable_chunks(
        queries,
        workdir=tmp_path / "serial",
        identity=identity,
        runner=_PickleableChunkRunner(reverse_completion=True),
        resume=True,
        force_recompute=False,
    )
    parallel = execute_restartable_chunks(
        queries,
        workdir=tmp_path / executor,
        identity=identity,
        runner=_PickleableChunkRunner(reverse_completion=True),
        resume=True,
        force_recompute=False,
        max_workers=2,
        executor=executor,
    )
    pd.testing.assert_frame_equal(serial.df, parallel.df)
    assert parallel.df["unique_sequence_id"].tolist() == queries["unique_sequence_id"].tolist()
    assert parallel.metadata["max_workers"] == 2
    assert parallel.metadata["executor"] == executor
    manifest = json.loads(Path(parallel.metadata["run_manifest_path"]).read_text())
    assert manifest["max_workers"] == 2


def test_parallel_cache_is_reusable_by_serial_execution(tmp_path: Path) -> None:
    identity = _identity(chunk_size=1)
    first = execute_restartable_chunks(
        _queries(),
        workdir=tmp_path,
        identity=identity,
        runner=_PickleableChunkRunner(),
        resume=True,
        force_recompute=False,
        max_workers=2,
        executor="thread",
    )
    second = execute_restartable_chunks(
        _queries(),
        workdir=tmp_path,
        identity=identity,
        runner=_PickleableChunkRunner(fail_index=0),
        resume=True,
        force_recompute=False,
        max_workers=1,
    )
    assert second.metadata["reused_chunk_count"] == first.metadata["chunk_count"]
    pd.testing.assert_frame_equal(first.df, second.df)


@pytest.mark.parametrize("max_workers", [0, -1, True, 1.5])
def test_parallel_worker_validation(tmp_path: Path, max_workers: Any) -> None:
    with pytest.raises(ValueError, match="positive integer"):
        execute_restartable_chunks(
            _queries(),
            workdir=tmp_path,
            identity=_identity(),
            runner=_PickleableChunkRunner(),
            resume=True,
            force_recompute=False,
            max_workers=max_workers,
        )


def test_parallel_failure_names_chunk_and_keeps_completed_manifests(tmp_path: Path) -> None:
    with pytest.raises(RuntimeError, match=r"index 1.*parallel injected failure"):
        execute_restartable_chunks(
            _queries(),
            workdir=tmp_path,
            identity=_identity(chunk_size=1),
            runner=_PickleableChunkRunner(fail_index=1),
            resume=True,
            force_recompute=False,
            max_workers=2,
            executor="thread",
        )
    manifests = sorted(tmp_path.glob("runs/*/chunks/*/manifest.json"))
    states = [json.loads(path.read_text())["status"] for path in manifests]
    assert "failed" in states
    assert all(state in {"complete", "failed", "running"} for state in states)
