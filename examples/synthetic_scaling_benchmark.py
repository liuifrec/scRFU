#!/usr/bin/env python3
"""Deterministic pure-Python scaling benchmark with no external RFU assets."""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import resource
import time
from collections.abc import Callable
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from scrfu import __version__, pp, tl

_AA = "ACDEFGHIKLMNPQRSTVWY"


def _aa_code(value: int, width: int = 8) -> str:
    letters: list[str] = []
    for _ in range(width):
        value, remainder = divmod(value, len(_AA))
        letters.append(_AA[remainder])
    return "C" + "".join(letters) + "F"


def generate_receptors(
    row_count: int,
    *,
    unique_fraction: float = 0.4,
    clone_distribution: str = "zipf",
    chain_mix: float = 0.9,
    missing_fraction: float = 0.02,
    sample_count: int = 100,
    donor_count: int = 20,
    phenotype_count: int = 8,
    random_state: int = 0,
) -> pd.DataFrame:
    """Generate a deterministic receptor table for stress and scaling tests."""
    if row_count < 1:
        raise ValueError("row_count must be positive.")
    if not 0 < unique_fraction <= 1:
        raise ValueError("unique_fraction must be in (0, 1].")
    if clone_distribution not in {"zipf", "uniform"}:
        raise ValueError("clone_distribution must be 'zipf' or 'uniform'.")
    if not 0 <= chain_mix <= 1 or not 0 <= missing_fraction < 1:
        raise ValueError("chain_mix and missing_fraction must be valid fractions.")
    rng = np.random.default_rng(random_state)
    unique_count = max(1, min(row_count, round(row_count * unique_fraction)))
    sequence_index = np.arange(unique_count, dtype=np.int64)
    if unique_count < row_count:
        remaining = row_count - unique_count
        probabilities = None
        if clone_distribution == "zipf":
            probabilities = 1.0 / np.arange(1, unique_count + 1, dtype=float)
            probabilities /= probabilities.sum()
        sequence_index = np.concatenate(
            [
                sequence_index,
                rng.choice(unique_count, size=remaining, replace=True, p=probabilities),
            ]
        )
        rng.shuffle(sequence_index)
    sequences = np.asarray([_aa_code(int(value)) for value in range(unique_count)], dtype=object)
    cdr3aa = sequences[sequence_index]
    chain = np.where(rng.random(row_count) < chain_mix, "TRB", "TRA").astype(object)
    v_index = sequence_index % 30 + 1
    v_call = np.asarray([f"TRBV{value}" for value in v_index], dtype=object)
    missing_v = rng.random(row_count) < missing_fraction
    v_call[missing_v] = pd.NA
    malformed = rng.random(row_count) < missing_fraction / 5
    cdr3aa = cdr3aa.astype(object, copy=True)
    cdr3aa[malformed] = "MALFORMED"
    sample_index = np.arange(row_count) % sample_count
    return pd.DataFrame(
        {
            "cell_id": [f"cell_{index:09d}" for index in range(row_count)],
            "chain": chain,
            "cdr3aa": cdr3aa,
            "v_call": v_call,
            "productive": rng.random(row_count) >= missing_fraction,
            "source_adapter": "synthetic_scaling",
            "source_row_id": [f"source_{index:09d}" for index in range(row_count)],
            "sample": [f"sample_{value:03d}" for value in sample_index],
            "donor": [f"donor_{value % donor_count:03d}" for value in sample_index],
            "phenotype": [f"state_{value % phenotype_count:02d}" for value in sequence_index],
        }
    )


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _peak_rss_bytes() -> int:
    factor = 1024 if platform.system() != "Darwin" else 1
    return int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * factor)


def _timed(operation: str, function: Callable[[], Any]) -> tuple[Any, dict[str, Any]]:
    wall_start = time.perf_counter()
    cpu_start = time.process_time()
    value = function()
    return value, {
        "operation": operation,
        "wall_seconds": time.perf_counter() - wall_start,
        "cpu_seconds": time.process_time() - cpu_start,
        "peak_rss_bytes": _peak_rss_bytes(),
    }


def benchmark_size(row_count: int, **generator_parameters: Any) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    receptors, timing = _timed(
        "generate", lambda: generate_receptors(row_count, **generator_parameters)
    )
    rows.append(timing)
    canonical, timing = _timed(
        "canonical_validation",
        lambda: pp.canonicalize_receptor_table(receptors, source_adapter="synthetic_scaling"),
    )
    pp.validate_receptor_table(canonical)
    rows.append(timing)
    unique, timing = _timed(
        "deduplication", lambda: canonical.drop_duplicates("cdr3aa", keep="first").copy()
    )
    rows.append(timing)
    unique["rfu_label"] = [f"RFU{index % 5000 + 1}" for index in range(len(unique))]
    unique["rfu_score"] = (np.arange(len(unique)) % 1000) / 1000
    unique["pass_thr"] = unique["rfu_score"].ge(0.2)
    assigned, timing = _timed(
        "reconstruction",
        lambda: canonical.merge(
            unique[["cdr3aa", "rfu_label", "rfu_score", "pass_thr"]],
            on="cdr3aa",
            how="left",
            sort=False,
            validate="many_to_one",
        ),
    )
    if (
        len(assigned) != row_count
        or assigned["input_row_id"].tolist() != canonical["input_row_id"].tolist()
    ):
        raise RuntimeError("Synthetic reconstruction changed input row count or order.")
    rows.append(timing)
    repertoire, timing = _timed(
        "repertoire_metrics",
        lambda: tl.repertoire_metrics(assigned, groupby="sample", weighting="cell"),
    )
    rows.append(
        timing | {"output_rows": len(repertoire), "output_columns": len(repertoire.columns)}
    )
    pseudobulk, timing = _timed(
        "pseudobulk",
        lambda: tl.rfu_pseudobulk(
            assigned, sample_key="sample", weighting="cell", normalize="count"
        ),
    )
    rows.append(
        timing
        | {"output_rows": pseudobulk.matrix.shape[0], "output_columns": pseudobulk.matrix.shape[1]}
    )
    overlap, timing = _timed("overlap", lambda: tl.rfu_overlap(pseudobulk, metric="cosine"))
    rows.append(
        timing | {"output_rows": len(overlap.matrix), "output_columns": len(overlap.matrix.columns)}
    )
    coupling, timing = _timed(
        "phenotype_coupling",
        lambda: tl.rfu_phenotype_coupling(
            assigned, phenotype_key="phenotype", sample_key="sample", weighting="cell"
        ),
    )
    rows.append(timing | {"output_rows": len(coupling), "output_columns": len(coupling.columns)})
    for row in rows:
        row.update(
            {
                "input_rows": row_count,
                "unique_sequences": int(canonical["cdr3aa"].nunique()),
                "deduplication_ratio": row_count / max(1, int(canonical["cdr3aa"].nunique())),
                **generator_parameters,
            }
        )
    return rows


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rows", type=int, nargs="+", default=[10_000, 100_000, 500_000])
    parser.add_argument("--unique-fraction", type=float, default=0.4)
    parser.add_argument("--clone-distribution", choices=["zipf", "uniform"], default="zipf")
    parser.add_argument("--chain-mix", type=float, default=0.9)
    parser.add_argument("--missing-fraction", type=float, default=0.02)
    parser.add_argument("--random-state", type=int, default=20260825)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args(argv)
    parameters = {
        "unique_fraction": args.unique_fraction,
        "clone_distribution": args.clone_distribution,
        "chain_mix": args.chain_mix,
        "missing_fraction": args.missing_fraction,
        "random_state": args.random_state,
    }
    results = pd.DataFrame(
        [row for size in args.rows for row in benchmark_size(size, **parameters)]
    )
    args.outdir.mkdir(parents=True, exist_ok=True)
    table_path = args.outdir / "synthetic_scaling.tsv"
    results.to_csv(table_path, sep="\t", index=False)
    manifest = {
        "schema_version": 1,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "scrfu_version": __version__,
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "parameters": parameters | {"row_counts": args.rows},
        "output": {
            "filename": table_path.name,
            "sha256": _sha256(table_path),
            "row_count": len(results),
        },
    }
    (args.outdir / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
