"""Acquire and audit inventoried HyPR-HN processed receptors, without RFU endpoints.

Only matched longitudinal TCR contigs and GEX barcodes are downloaded (explicit
--download). Matrices and other cohorts are excluded. Missing files stay unknown.
This stage freezes the development configuration and establishes actual coverage
before any external RFU assignment or outcome inspection.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlparse

import pandas as pd

from scrfu.adapters import adapt_cellranger_vdj
from scrfu.io import file_sha256

from .radiation_methods import REPO, completed, json_save, save, validate_boundary


def select_files(inventory: pd.DataFrame, pairs: pd.DataFrame) -> pd.DataFrame:
    keys = ["donor", "visit", "compartment"]
    approved = pairs[pairs.paired_processed_GEX_TCR_available.eq(True)][keys]
    if approved.duplicated(keys).any():
        raise ValueError("Ambiguous donor/visit/compartment pairing.")
    selected = inventory[inventory.longitudinal_hypr_sample.eq(True)].merge(
        approved, on=keys, validate="many_to_one"
    )
    selected = selected[
        selected.url.str.endswith(("filtered_contig_annotations.csv.gz", "barcodes.tsv.gz"))
    ].copy()
    if selected.duplicated(keys + ["assay"]).any():
        raise ValueError("Multiple processed inputs for one sample/assay.")
    for _, group in selected.groupby(keys):
        if set(group.assay) != {"GEX", "TCR"}:
            raise ValueError("Processed pairing lacks an assay.")
    for url in selected.url:
        parts = urlparse(url)
        if parts.scheme != "https" or parts.hostname != "ftp.ncbi.nlm.nih.gov":
            raise ValueError("Only the inventoried HTTPS NCBI files are accepted.")
    return selected.sort_values(keys + ["assay"]).reset_index(drop=True)


def validate_gzip(path: Path) -> None:
    with gzip.open(path, "rb") as stream:
        while stream.read(1024 * 1024):
            pass


def obtain(url: str, directory: Path, download: bool) -> dict:
    path = directory / Path(urlparse(url).path).name
    result = {"url": url, "path": str(path.resolve())}
    try:
        if not path.exists():
            if not download:
                return {**result, "status": "missing_not_downloaded"}
            partial = path.with_name(path.name + ".part")
            subprocess.run(
                [
                    "curl",
                    "-fL",
                    "--silent",
                    "--show-error",
                    "--connect-timeout",
                    "20",
                    "--max-time",
                    "150",
                    "--retry",
                    "2",
                    "--retry-delay",
                    "3",
                    "--retry-max-time",
                    "180",
                    "--continue-at",
                    "-",
                    "--output",
                    str(partial),
                    url,
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            validate_gzip(partial)
            partial.rename(path)
        validate_gzip(path)
        return {
            **result,
            "status": "gzip_verified",
            "bytes": path.stat().st_size,
            "sha256": file_sha256(path),
            "checksum_scope": "locally_computed_SHA256_and_gzip_CRC; no_publisher_hash_supplied_in_inventory",
        }
    except (OSError, EOFError, subprocess.CalledProcessError) as error:
        return {**result, "status": "access_or_integrity_failure", "error": str(error)[:300]}


def treatment_label(characteristics: str, visit: int) -> str:
    match = re.search(r"(?:^|;\s*)treatment:\s*(.+?)(?:;|$)", characteristics)
    label = match.group(1).strip() if match else ""
    # Preserve the documented GEO typo in provenance; normalize only for validation.
    normalized = label.replace("Radition", "Radiation")
    expected = {1: "Pre-Treatment", 2: "Last Day of Radiation", 3: "6 Weeks Post-Radiation"}
    if normalized != expected.get(visit):
        raise ValueError("Visit number disagrees with source treatment metadata.")
    return normalized


def pair_receptors(
    contigs: pd.DataFrame, barcodes: pd.Series, sample_id: str
) -> tuple[pd.DataFrame, dict]:
    if barcodes.isna().any() or barcodes.duplicated().any():
        raise ValueError("GEX barcode list must be unique and nonmissing.")
    all_trb = adapt_cellranger_vdj(contigs, chain="TRB", primary_chain=False).receptors
    primary = adapt_cellranger_vdj(contigs, chain="TRB", primary_chain=True).receptors
    selected = primary[primary.cell_id.isin(barcodes)].copy()
    qc = {
        "gex_barcodes": len(barcodes),
        "contig_rows": len(contigs),
        "filtered_productive_TRB_contigs": len(all_trb),
        "primary_TRB_cells": len(primary),
        "primary_TRB_cells_in_GEX": len(selected),
        "primary_TRB_cells_not_in_GEX": len(primary) - len(selected),
        "cells_with_multiple_productive_TRB": int(all_trb.groupby("cell_id").size().gt(1).sum()),
        "primary_chain_policy": "existing_adapter: productive/high_confidence; highest_UMI_then_reads_then_source_order",
        "RNA_state_labels_available": False,
    }
    selected["source_barcode"] = selected.cell_id
    selected["cell_id"] = sample_id + ":" + selected.cell_id
    selected["input_row_id"] = sample_id + ":" + selected.input_row_id
    selected["sample_id"] = sample_id
    return selected, qc


def run(args: argparse.Namespace) -> None:
    if args.out.resolve().is_relative_to(REPO) or args.raw.resolve().is_relative_to(REPO):
        raise ValueError("External biological inputs and outputs must stay outside the repository.")
    config = json.loads(args.config.read_text())
    inventory, pairs = (pd.read_csv(path, sep="\t") for path in (args.inventory, args.pairs))
    selected = select_files(inventory, pairs)
    args.raw.mkdir(parents=True, exist_ok=True)
    args.out.mkdir(parents=True, exist_ok=True)
    # Configuration snapshot precedes any downloaded receptor inspection.
    config_path = args.out / "frozen_development_configuration.json"
    if config_path.exists() and json.loads(config_path.read_text()) != config:
        raise ValueError("Frozen external configuration changed.")
    json_save(config, config_path)
    with ThreadPoolExecutor(max_workers=2) as pool:
        downloads = list(pool.map(lambda u: obtain(u, args.raw, args.download), selected.url))
    files = pd.DataFrame(downloads)
    save(files, args.out, "processed_file_manifest.tsv")
    identity = {
        "inventory_sha256": file_sha256(args.inventory),
        "pairs_sha256": file_sha256(args.pairs),
        "config_sha256": file_sha256(args.config),
        "script_sha256": file_sha256(Path(__file__)),
        "adapter_sha256": file_sha256(REPO / "src/scrfu/adapters.py"),
        "inputs": {x["path"]: x.get("sha256") for x in downloads},
    }
    fingerprint = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
    if completed(args.out, fingerprint):
        print("GSE280982 prepared-source checkpoint verified and reused.")
        return
    json_save({"status": "running", "fingerprint": fingerprint}, args.out / "completion.json")
    rows, receptor_tables = [], []
    files_by_url = files.set_index("url").to_dict(orient="index")
    for row in pairs.to_dict(orient="records"):
        donor, visit, compartment = int(row["donor"]), int(row["visit"]), row["compartment"]
        sample = f"HP{donor:02}_{compartment}_V{visit}"
        meta = {
            "sample_id": sample,
            "donor": f"HP{donor:02}",
            "visit": visit,
            "compartment": compartment,
            "dataset_id": "GSE280982",
            "authorized": True,
        }
        if not row["paired_processed_GEX_TCR_available"]:
            rows.append({**meta, "status": "no_released_paired_TCR"})
            continue
        source = selected[
            selected.donor.eq(donor)
            & selected.visit.eq(visit)
            & selected.compartment.eq(compartment)
        ].set_index("assay")
        treatments = [treatment_label(str(c), visit) for c in source.characteristics]
        if len(set(treatments)) != 1:
            raise ValueError("TCR and GEX treatment annotations disagree.")
        meta["treatment_time"] = treatments[0]
        if any(files_by_url[u]["status"] != "gzip_verified" for u in source.url):
            rows.append({**meta, "status": "unknown_due_to_unavailable_file"})
            continue
        contigs = pd.read_csv(files_by_url[source.loc["TCR", "url"]]["path"])
        barcodes = pd.read_csv(
            files_by_url[source.loc["GEX", "url"]]["path"], sep="\t", header=None
        )[0]
        receptors, qc = pair_receptors(contigs, barcodes, sample)
        receptors = receptors.assign(**meta)
        receptor_tables.append(receptors)
        rows.append(
            {
                **meta,
                **qc,
                "status": "prepared_before_RFU_assignment",
                "TCR_accession": source.loc["TCR", "accession"],
                "GEX_accession": source.loc["GEX", "accession"],
                "minimum_cells_before_RFU_filter": len(receptors) >= config["min_cells_per_visit"],
            }
        )
    registry = pd.DataFrame(rows)
    save(registry, args.out, "gse280982_verified_sample_registry.tsv")
    if receptor_tables:
        receptors = pd.concat(receptor_tables, ignore_index=True)
        validate_boundary(receptors, registry)
        if receptors.cell_id.duplicated().any():
            raise ValueError("Cell identifiers collide across samples.")
        receptors.to_parquet(args.out / "gse280982_primary_trb_matched_gex.parquet", index=False)
    else:
        receptors = pd.DataFrame()
    summary = {
        "downloaded_files": int(files.status.eq("gzip_verified").sum()),
        "required_files": len(files),
        "bytes": int(files.bytes.fillna(0).sum()) if "bytes" in files else 0,
        "available_tumor_visits": int(
            registry.compartment.eq("tumor")
            .mul(registry.status.eq("prepared_before_RFU_assignment"))
            .sum()
        ),
        "available_blood_visits": int(
            registry.compartment.eq("blood")
            .mul(registry.status.eq("prepared_before_RFU_assignment"))
            .sum()
        ),
        "matched_primary_TRB_cells": len(receptors),
        "external_RFUs_assigned": False,
        "external_endpoints_computed": False,
        "cell_state_blocker": "GEX matrices/barcodes provide no author cell-state labels; no new RNA clustering undertaken",
    }
    json_save(summary, args.out / "evidence_counts.json")
    json_save(
        {
            **identity,
            "command": sys.argv,
            "python": sys.executable,
            "git_sha": subprocess.check_output(
                ["git", "rev-parse", "HEAD"], cwd=REPO, text=True
            ).strip(),
            "prior_access": "metadata/publication audited before external endpoints; not an untouched dataset",
            "timing_source": "https://doi.org/10.1038/s41467-025-60827-w and pinned GEO sample characteristics",
            "unknown_vs_zero": "Unavailable TCR entries and files have missing counts, never zeros",
        },
        args.out / "provenance.json",
    )
    status = "complete" if files.status.eq("gzip_verified").all() else "blocked_acquisition"
    json_save(
        {
            "status": status,
            "fingerprint": fingerprint,
            "outputs": {
                p.name: file_sha256(p)
                for p in args.out.iterdir()
                if p.is_file() and p.name != "completion.json"
            },
        },
        args.out / "completion.json",
    )
    print(json.dumps(summary, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ("inventory", "pairs", "raw", "out"):
        parser.add_argument("--" + key, type=Path, required=True)
    parser.add_argument(
        "--config", type=Path, default=REPO / "manuscript/config/radiation_methods_v1.json"
    )
    parser.add_argument("--download", action="store_true")
    run(parser.parse_args())


if __name__ == "__main__":
    main()
