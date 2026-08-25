#!/usr/bin/env python3
"""Run the public synthetic scRFU tutorial with a mock or external RFU backend."""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
from datetime import datetime, timezone
from importlib.resources import files
from pathlib import Path
from typing import Any

import pandas as pd

from scrfu import __version__, bcr, pp, tl


def _fixture(name: str) -> Path:
    return Path(str(files("scrfu").joinpath("data", name)))


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _mock_assignments(receptors: pd.DataFrame) -> pd.DataFrame:
    frame = receptors.copy()
    codes, _ = pd.factorize(frame["cdr3aa"], sort=False)
    frame["unique_sequence_id"] = [f"tutorial_sequence_{value:04d}" for value in codes]
    frame["rfu_id"] = frame.pop("mock_rfu_id")
    frame["rfu_label"] = frame.pop("mock_rfu_label")
    frame["rfu_score"] = frame.pop("mock_rfu_score")
    frame["pass_thr"] = frame.pop("mock_pass_thr").astype("boolean")
    frame["eligibility_status"] = "eligible"
    return frame


def run_tutorial(
    *,
    outdir: Path,
    backend: str = "mock",
    rfu_dir: Path | None = None,
) -> dict[str, Any]:
    receptor_path = _fixture("tutorial_receptors.tsv")
    vdjdb_path = _fixture("tutorial_vdjdb.tsv")
    bcr_path = _fixture("tutorial_bcr.tsv")
    receptors = pp.canonicalize_receptor_table(pd.read_csv(receptor_path, sep="\t"))
    receptor_qc = pp.validate_receptor_table(receptors)
    if backend == "mock":
        assigned = _mock_assignments(receptors)
        backend_notice = (
            "Synthetic assignments for software demonstration; not official RFU results."
        )
        rfu_provenance: dict[str, Any] = {"backend": "synthetic_mock", "external_rfu": False}
    elif backend == "rfu_repo":
        if rfu_dir is None:
            raise ValueError("--rfu-dir is required with --backend rfu_repo.")
        result = tl.call_rfu_table(
            receptors.drop(columns=[column for column in receptors if column.startswith("mock_")]),
            rfu_dir=rfu_dir,
            chunk_size=4,
            workdir=outdir / "rfu_backend",
        )
        assigned = result.per_row.merge(
            receptors[["input_row_id", "sample", "donor", "time", "phenotype", "clonotype_id"]],
            on="input_row_id",
            how="left",
            sort=False,
            validate="one_to_one",
        )
        backend_notice = "Official external RFU checkout requested by the user."
        rfu_provenance = result.provenance
    else:
        raise ValueError("backend must be 'mock' or 'rfu_repo'.")

    repertoire = tl.repertoire_metrics(assigned, groupby="sample", weighting="cell", chain="TRB")
    pseudobulk = tl.rfu_pseudobulk(
        assigned,
        sample_key="sample",
        phenotype_keys=["donor", "time"],
        weighting="cell",
        normalize="proportion",
    )
    coupling = tl.rfu_phenotype_coupling(
        assigned, phenotype_key="phenotype", sample_key="sample", weighting="cell"
    )
    longitudinal = tl.rfu_longitudinal_matrix(
        assigned,
        sample_key="sample",
        donor_key="donor",
        time_key="time",
        weighting="cell",
        normalize="proportion",
    )
    similarity = tl.longitudinal_similarity(longitudinal, metric="cosine")
    retrieval = tl.donor_retrieval(
        longitudinal, metric="cosine", top_k=3, exclude_same_timepoint=True
    )
    reference = tl.load_vdjdb_reference(vdjdb_path, release_label="synthetic-tutorial-v1")
    evidence = tl.annotate_vdjdb(assigned, reference, match_mode="cdr3_v", chain="TRB")
    evidence_summary = tl.summarize_vdjdb_evidence(assigned, evidence)
    bcr_result = bcr.prepare_bcr_table(
        pd.read_csv(bcr_path, sep="\t"), source_label="synthetic_tutorial"
    )
    bcr_states = bcr.bcr_state_features(bcr_result.receptors, bcr_result.pairs)
    bcr_features = bcr.bcr_feature_matrix(bcr_result.receptors, pairs=bcr_result.pairs)

    outdir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "assigned_receptors.tsv": assigned,
        "repertoire_metrics.tsv": repertoire,
        "rfu_pseudobulk_proportions.tsv": pseudobulk.matrix.reset_index(),
        "phenotype_coupling.tsv": coupling,
        "longitudinal_similarity.tsv": similarity,
        "donor_retrieval.tsv": retrieval,
        "synthetic_vdjdb_evidence.tsv": evidence,
        "synthetic_vdjdb_sequence_summary.tsv": evidence_summary.sequence_summary,
        "synthetic_vdjdb_row_summary.tsv": evidence_summary.row_summary,
        "bcr_receptors.tsv": bcr_result.receptors,
        "bcr_pairs.tsv": bcr_result.pairs,
        "bcr_state_features.tsv": bcr_states,
        "bcr_feature_matrix.tsv": bcr_features.features,
        "bcr_feature_missingness.tsv": bcr_features.missingness,
    }
    output_manifest: dict[str, Any] = {}
    for name, table in outputs.items():
        path = outdir / name
        table.to_csv(path, sep="\t", index=False)
        output_manifest[name] = {
            "sha256": _sha256(path),
            "rows": len(table),
            "columns": len(table.columns),
        }
    manifest = {
        "schema_version": 1,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "scrfu_version": __version__,
        "python_version": platform.python_version(),
        "backend": backend,
        "backend_notice": backend_notice,
        "rfu_provenance": rfu_provenance,
        "receptor_qc": receptor_qc,
        "bcr_qc": bcr_result.qc,
        "bcr_feature_parameters": bcr_features.parameters,
        "fixtures": {
            path.name: {"sha256": _sha256(path), "synthetic": True}
            for path in (receptor_path, vdjdb_path, bcr_path)
        },
        "outputs": output_manifest,
    }
    (outdir / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return manifest


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--backend", choices=["mock", "rfu_repo"], default="mock")
    parser.add_argument("--rfu-dir", type=Path)
    args = parser.parse_args(argv)
    manifest = run_tutorial(
        outdir=args.outdir.expanduser().resolve(),
        backend=args.backend,
        rfu_dir=args.rfu_dir,
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
