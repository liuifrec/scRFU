#!/usr/bin/env python3
"""Validate VDJdb linkage for the public Scirpy wu2020_3k native AIRR run.

The public dataset's cell identifiers encode the sample before the first
underscore; this convention is used only for the descriptive sample summary.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import awkward as ak
import mudata
import pandas as pd

import scrfu


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _flatten_native(airr: Any) -> tuple[pd.DataFrame, pd.DataFrame]:
    airr_records = ak.to_list(airr.obsm["airr"])
    rfu_records = ak.to_list(airr.obsm["scrfu"])
    if len(airr_records) != len(rfu_records):
        raise RuntimeError("AIRR and scRFU observation counts differ.")

    all_chains: list[dict[str, Any]] = []
    queries: list[dict[str, Any]] = []
    for obs_index, (cell_id, source_chains, result_chains) in enumerate(
        zip(airr.obs_names.astype(str), airr_records, rfu_records, strict=True)
    ):
        if len(source_chains) != len(result_chains):
            raise RuntimeError(f"AIRR/scRFU chain alignment differs for {cell_id!r}.")
        for chain_index, (source, result) in enumerate(
            zip(source_chains, result_chains, strict=True)
        ):
            row_id = f"native_{obs_index:09d}_{chain_index:03d}"
            chain_row = {
                "input_row_id": row_id,
                "cell_id": cell_id,
                "sample_id": cell_id.split("_", 1)[0],
                "chain_index": chain_index,
                "chain": source.get("locus"),
                "cdr3aa": source.get("junction_aa"),
                "v_call": source.get("v_call"),
                "eligible": bool(result.get("eligible", False)),
                "unique_sequence_id": result.get("unique_sequence_id"),
                "rfu_id": result.get("rfu_id"),
                "rfu_label": result.get("rfu_label"),
                "rfu_score": result.get("rfu_score"),
                "pass_thr": result.get("threshold_pass"),
                "assignment_status": result.get("assignment_status"),
            }
            all_chains.append(chain_row)
            if chain_row["eligible"]:
                queries.append(chain_row)
    return pd.DataFrame(all_chains), pd.DataFrame(queries)


def run(args: argparse.Namespace) -> dict[str, Any]:
    output = args.output_dir.expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    mdata = mudata.read_h5mu(args.input)
    validation = scrfu.tl.validate_scrfu_schema(mdata)
    airr = mdata.mod[args.airr_mod]
    all_chains, queries = _flatten_native(airr)
    if queries["input_row_id"].duplicated().any():
        raise RuntimeError("Native match queries contain duplicate row identifiers.")
    if queries["unique_sequence_id"].isna().any():
        raise RuntimeError("Eligible native chains are missing RFU sequence identities.")
    if queries.groupby("unique_sequence_id")["cdr3aa"].nunique().gt(1).any():
        raise RuntimeError("RFU sequence identity maps to conflicting CDR3 values.")

    reference = scrfu.tl.load_vdjdb_reference(
        args.vdjdb,
        release_label=args.release,
        expected_sha256=args.vdjdb_sha256,
        sep="\t",
    )
    output_tables: dict[str, pd.DataFrame] = {
        "native_chain_queries.tsv.gz": queries,
    }
    mode_metrics: dict[str, Any] = {}
    chain_summary = all_chains.loc[
        :, ["input_row_id", "cell_id", "sample_id", "chain_index", "chain", "eligible"]
    ].copy()
    for mode in ("cdr3", "cdr3_v"):
        evidence = scrfu.tl.annotate_vdjdb(
            queries,
            reference,
            match_mode=mode,
            chain="TRB",
            expand_rows=False,
        )
        summary = scrfu.tl.summarize_vdjdb_evidence(queries, evidence)
        if len(summary.row_summary) != len(queries):
            raise RuntimeError(f"{mode} row reconstruction changed the query count.")
        if summary.row_summary["input_row_id"].tolist() != queries["input_row_id"].tolist():
            raise RuntimeError(f"{mode} row reconstruction changed input order.")
        row_columns = [
            "input_row_id",
            "has_vdjdb_evidence",
            "evidence_record_count",
            "distinct_epitope_count",
            "antigen_ambiguity",
        ]
        row_summary = summary.row_summary.loc[:, row_columns].copy()
        row_summary = row_summary.add_prefix(f"{mode}_").rename(
            columns={f"{mode}_input_row_id": "input_row_id"}
        )
        before = len(chain_summary)
        chain_summary = chain_summary.merge(
            row_summary, on="input_row_id", how="left", sort=False, validate="one_to_one"
        )
        if len(chain_summary) != before:
            raise RuntimeError(f"{mode} evidence expanded the native chain table.")
        evidence_name = f"{mode}_vdjdb_matches_long.tsv.gz"
        row_name = f"{mode}_row_summary.tsv.gz"
        output_tables[evidence_name] = evidence
        output_tables[row_name] = summary.row_summary
        mode_metrics[mode] = {
            "evidence_rows": len(evidence),
            "matched_rfu_sequence_identities": int(evidence["unique_sequence_id"].nunique()),
            "matched_query_variants": int(evidence["match_query_id"].nunique()),
            "matched_input_chains": int(summary.row_summary["has_vdjdb_evidence"].sum()),
            "ambiguous_input_chains": int(summary.row_summary["antigen_ambiguity"].sum()),
            "row_reconstruction_count": len(summary.row_summary),
        }

    evidence_columns = [f"{mode}_has_vdjdb_evidence" for mode in ("cdr3", "cdr3_v")]
    chain_summary[evidence_columns] = (
        chain_summary[evidence_columns].astype("boolean").fillna(False).astype(bool)
    )
    if len(chain_summary) != len(all_chains):
        raise RuntimeError("Final VDJdb linkage changed the AIRR chain count.")
    cell_summary = (
        chain_summary.groupby(["cell_id", "sample_id"], sort=False, observed=True)
        .agg(
            airr_chain_count=("input_row_id", "size"),
            eligible_trb_count=("eligible", "sum"),
            cdr3_matched_chain_count=("cdr3_has_vdjdb_evidence", "sum"),
            cdr3_v_matched_chain_count=("cdr3_v_has_vdjdb_evidence", "sum"),
        )
        .reset_index()
    )
    cell_summary = pd.DataFrame({"cell_id": airr.obs_names.astype(str)}).merge(
        cell_summary, on="cell_id", how="left", sort=False, validate="one_to_one"
    )
    count_columns = [column for column in cell_summary if column.endswith("_count")]
    cell_summary[count_columns] = cell_summary[count_columns].fillna(0).astype(int)
    cell_summary["sample_id"] = cell_summary["cell_id"].str.split("_", n=1).str[0]
    sample_summary = (
        cell_summary.groupby("sample_id", sort=True, observed=True)
        .agg(
            cell_count=("cell_id", "size"),
            eligible_trb_count=("eligible_trb_count", "sum"),
            cdr3_matched_chain_count=("cdr3_matched_chain_count", "sum"),
            cdr3_v_matched_chain_count=("cdr3_v_matched_chain_count", "sum"),
        )
        .reset_index()
    )
    output_tables.update(
        {
            "native_chain_vdjdb_summary.tsv.gz": chain_summary,
            "native_cell_vdjdb_summary.tsv.gz": cell_summary,
            "native_sample_vdjdb_summary.tsv": sample_summary,
        }
    )
    for name, table in output_tables.items():
        compression = {"method": "gzip", "mtime": 0} if name.endswith(".gz") else None
        table.to_csv(output / name, sep="\t", index=False, compression=compression)

    manifest: dict[str, Any] = {
        "schema_version": 1,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "dataset": "Scirpy wu2020_3k native AIRR bounded VDJdb linkage",
        "public_source": "https://scverse-exampledata.s3.eu-west-1.amazonaws.com/scirpy/wu2020_3k.h5mu",
        "scrfu_version": scrfu.__version__,
        "python_version": platform.python_version(),
        "parameters": {
            "match_modes": ["cdr3", "cdr3_v"],
            "chain": "TRB",
            "v_gene_normalization": "strip_allele",
            "evidence_representation": "unexpanded authoritative long table",
            "rfu_identity": "canonical exact CDR3",
        },
        "inputs": {
            "native_h5mu_sha256": _sha256(args.input),
            "native_h5mu_size_bytes": args.input.stat().st_size,
            "vdjdb_release": args.release,
            "vdjdb_sha256": reference.provenance["sha256"],
            "vdjdb_rows": len(reference.table),
            "rfu_reference_identity": validation["reference_identity"],
        },
        "counts": {
            "cells": airr.n_obs,
            "airr_chains": len(all_chains),
            "eligible_trb_chains": len(queries),
            "rfu_sequence_identities": int(queries["unique_sequence_id"].nunique()),
            **mode_metrics,
        },
        "invariants": {
            "airr_chain_count_preserved": len(chain_summary) == len(all_chains),
            "cell_count_preserved": len(cell_summary) == airr.n_obs,
            "rfu_identity_cdr3_consistent": True,
            "row_order_preserved": True,
        },
        "outputs": {
            name: {
                "sha256": _sha256(output / name),
                "size_bytes": (output / name).stat().st_size,
                "rows": len(table),
            }
            for name, table in output_tables.items()
        },
    }
    (output / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Native H5MU containing aligned obsm['airr'] and obsm['scrfu'] records.",
    )
    parser.add_argument(
        "--vdjdb", type=Path, required=True, help="Pinned external VDJdb TSV reference."
    )
    parser.add_argument("--release", required=True, help="Explicit VDJdb release label.")
    parser.add_argument(
        "--vdjdb-sha256",
        required=True,
        help="Expected SHA256 of the pinned VDJdb TSV reference.",
    )
    parser.add_argument(
        "--airr-mod", default="airr", help="MuData modality containing native AIRR records."
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Explicit external directory for linkage tables and the run manifest.",
    )
    args = parser.parse_args()
    print(json.dumps(run(args), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
