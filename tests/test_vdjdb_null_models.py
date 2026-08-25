from __future__ import annotations

import json

import pandas as pd

from examples.vdjdb_null_models import run_null_models
from scrfu.tl import annotate_vdjdb, load_vdjdb_reference


def test_null_model_workflow_writes_four_prespecified_models(tmp_path) -> None:
    sequences = pd.DataFrame(
        {
            "input_row_id": ["a", "b", "c", "d"],
            "unique_sequence_id": ["a", "b", "c", "d"],
            "cdr3aa": ["CASSAAA", "CASSAAD", "CASSCCC", "CASSCCD"],
            "v_call": ["TRBV1", "TRBV1", "TRBV2", "TRBV2"],
            "chain": ["TRB"] * 4,
            "rfu_label": ["R1", "R1", "R2", "R2"],
            "pass_thr": [True] * 4,
        }
    )
    reference = load_vdjdb_reference(
        sequences[["cdr3aa", "v_call", "chain"]].assign(epitope=["A", "A", "B", "B"]),
        release_label="synthetic",
    )
    evidence = annotate_vdjdb(sequences, reference, expand_rows=False)
    sequence_path = tmp_path / "sequences.tsv"
    evidence_path = tmp_path / "evidence.tsv"
    outdir = tmp_path / "nulls"
    sequences.to_csv(sequence_path, sep="\t", index=False)
    evidence.to_csv(evidence_path, sep="\t", index=False)

    manifest = run_null_models(
        rfu_sequences=sequence_path,
        evidence=evidence_path,
        outdir=outdir,
        assignment_policy="nearest",
        ambiguity_policy="fractional",
        n_permutations=5,
        random_state=4,
        save_values=True,
    )

    summary = pd.read_csv(outdir / "antigen_null_model_summary.tsv", sep="\t")
    assert summary["null_model"].tolist() == [
        "unrestricted",
        "cdr3_length",
        "trbv",
        "trbv_cdr3_length",
    ]
    assert summary["n_permutations"].eq(5).all()
    assert len(pd.read_csv(outdir / "antigen_null_model_values.tsv.gz", sep="\t")) == 20
    saved = json.loads((outdir / "null_model_manifest.json").read_text())
    assert saved["parameters"]["aggregation"] == "unique_sequence_antigen"
    assert manifest["inputs"]["rfu_sequence_rows"] == 4
