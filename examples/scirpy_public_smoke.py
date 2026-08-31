#!/usr/bin/env python3
"""Run a small public Scirpy AIRR interoperability smoke with an external RFU checkout."""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import time
from datetime import datetime, timezone
from pathlib import Path

import awkward as ak
import mudata
import pandas as pd
import scirpy

import scrfu


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run(input_path: Path, output_dir: Path, rfu_dir: Path) -> dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=True)
    mdata = mudata.read_h5mu(input_path)
    if "airr" not in mdata.mod:
        raise KeyError("Public smoke input has no 'airr' modality.")
    airr = mdata.mod["airr"]
    before = ak.to_list(airr.obsm["airr"])
    started = time.perf_counter()
    explicit = scrfu.tl.assign_rfu(
        mdata,
        rfu_dir=rfu_dir,
        chunk_size=5000,
        max_workers=1,
        workdir=output_dir / "rfu_cache",
        cell_summary=True,
        summary_policy="ambiguity_aware",
        inplace=False,
    )
    if explicit is None:
        raise RuntimeError("Explicit native assignment returned no structured result.")
    native_records = ak.to_list(explicit.chain_records)
    parity_mismatches = 0
    for row in explicit.table_result.per_row.to_dict(orient="records"):
        _, obs_text, chain_text = str(row["input_row_id"]).split("_")
        chain_record = native_records[int(obs_text)][int(chain_text)]
        values = (
            chain_record["rfu_id"],
            chain_record["rfu_label"],
            chain_record["rfu_score"],
            chain_record["threshold_pass"],
        )
        expected = (row["rfu_id"], row["rfu_label"], row["rfu_score"], row["pass_thr"])
        if not all(
            (pd.isna(left) and pd.isna(right)) or left == right
            for left, right in zip(values, expected, strict=True)
        ):
            parity_mismatches += 1
    scrfu.tl.assign_rfu(
        mdata,
        rfu_dir=rfu_dir,
        chunk_size=5000,
        max_workers=1,
        workdir=output_dir / "rfu_cache",
        cell_summary=True,
        summary_policy="ambiguity_aware",
    )
    elapsed = time.perf_counter() - started
    if ak.to_list(airr.obsm["airr"]) != before:
        raise RuntimeError("Native assignment mutated public AIRR records.")
    validation = scrfu.tl.validate_scrfu_schema(mdata)
    annotated = output_dir / "wu2020_3k_scrfu.h5mu"
    mdata.write_h5mu(annotated)
    reloaded = mudata.read_h5mu(annotated)
    reload_validation = scrfu.tl.validate_scrfu_schema(reloaded)
    chain_records = ak.to_list(reloaded.mod["airr"].obsm["scrfu"])
    trb = [chain for cell in chain_records for chain in cell if chain["locus"] == "TRB"]
    manifest: dict[str, object] = {
        "schema_version": 1,
        "dataset": "Scirpy wu2020_3k",
        "public_source": "https://scverse-exampledata.s3.eu-west-1.amazonaws.com/scirpy/wu2020_3k.h5mu",
        "input_sha256": _sha256(input_path),
        "input_size_bytes": input_path.stat().st_size,
        "output_sha256": _sha256(annotated),
        "output_size_bytes": annotated.stat().st_size,
        "scrfu_version": scrfu.__version__,
        "scirpy_version": scirpy.__version__,
        "python_version": platform.python_version(),
        "created_at": datetime.now(timezone.utc).isoformat(),
        "elapsed_seconds": elapsed,
        "cells": airr.n_obs,
        "airr_chains": validation["chain_count"],
        "trb_chains": len(trb),
        "eligible_trb_chains": sum(bool(chain["eligible"]) for chain in trb),
        "threshold_qualified_trb_chains": sum(bool(chain["threshold_pass"]) for chain in trb),
        "canonical_table_parity_rows": len(explicit.table_result.per_row),
        "canonical_table_parity_mismatches": parity_mismatches,
        "reference_identity": validation["reference_identity"],
        "reload_validation": reload_validation,
        "assignment_parameters": {
            "mode": "standard",
            "threshold": 0.6,
            "deduplicate": True,
            "chunk_size": 5000,
            "max_workers": 1,
            "summary_policy": "ambiguity_aware",
        },
    }
    path = output_dir / "run_manifest.json"
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    pd.DataFrame(
        {
            "metric": ["cells", "airr_chains", "trb_chains", "eligible_trb_chains"],
            "value": [
                airr.n_obs,
                validation["chain_count"],
                len(trb),
                manifest["eligible_trb_chains"],
            ],
        }
    ).to_csv(output_dir / "summary.tsv", sep="\t", index=False)
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--rfu-dir", type=Path, required=True)
    args = parser.parse_args()
    print(
        json.dumps(
            run(args.input.resolve(), args.output_dir.resolve(), args.rfu_dir.resolve()),
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
