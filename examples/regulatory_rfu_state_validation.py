"""Donor-held-out RFU versus receptor-gene prediction using cached GSE190905.

No assignment or outcome-derived cell labeling occurs here. See the frozen
specification in docs/regulatory_execution_state.md. Biological outputs are local.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, log_loss
from sklearn.preprocessing import OneHotEncoder


def sha256(path: Path) -> str:
    hasher = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def split_donor(
    frame: pd.DataFrame, donor: str, *, purge_shared: bool = True
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Donor generalization, optionally restricted to unseen receptor identities."""
    test = frame.loc[frame.patient == donor].copy()
    keep = frame.patient != donor
    if purge_shared:
        keep &= ~frame.clone_key.isin(test.clone_key)
    train = frame.loc[keep].copy()
    assert not set(train.patient) & set(test.patient)
    if purge_shared:
        assert not set(train.clone_key) & set(test.clone_key)
    return train, test


def paired_deltas(metrics: pd.DataFrame) -> pd.DataFrame:
    """Extended minus baseline: positive log-loss change means deterioration."""
    paired = metrics.pivot(
        index=["held_out_donor", "weighting"], columns="model", values="log_loss"
    )
    if paired.isna().any().any():
        raise ValueError("Unpaired model evaluation")
    result = paired.rename(
        columns={
            "TRBV_TRBJ_length": "baseline_log_loss",
            "TRBV_TRBJ_length_RFU": "extended_log_loss",
        }
    ).reset_index()
    result["delta_log_loss"] = result.extended_log_loss - result.baseline_log_loss
    return result


def audit_predictions(predictions: pd.DataFrame, metrics: pd.DataFrame) -> None:
    """Recalculate scores and verify identical evaluation cells, labels and weights."""
    for (donor, mode), subset in predictions.groupby(["held_out_donor", "weighting"]):
        groups = {
            name: sub.sort_values("cell_id").set_index("cell_id")
            for name, sub in subset.groupby("model")
        }
        if set(groups) != {"TRBV_TRBJ_length", "TRBV_TRBJ_length_RFU"}:
            raise ValueError("Missing paired predictions")
        base, extended = groups["TRBV_TRBJ_length"], groups["TRBV_TRBJ_length_RFU"]
        if not base.index.is_unique or not extended.index.is_unique:
            raise ValueError("Duplicated evaluation observations")
        pd.testing.assert_frame_equal(
            base[["cluster", "evaluation_weight"]], extended[["cluster", "evaluation_weight"]]
        )
        for model_name, group in groups.items():
            columns = sorted(c for c in group if c.startswith("prob:") and group[c].notna().any())
            probabilities = group[columns].to_numpy()
            if not np.allclose(probabilities.sum(axis=1), 1):
                raise ValueError("Probabilities do not sum to one")
            labels = [c.removeprefix("prob:") for c in columns]
            score = log_loss(
                group.cluster, probabilities, labels=labels, sample_weight=group.evaluation_weight
            )
            saved = metrics.loc[
                (metrics.held_out_donor == donor)
                & (metrics.weighting == mode)
                & (metrics.model == model_name),
                "log_loss",
            ]
            if len(saved) != 1 or not np.isclose(score, saved.iloc[0], rtol=0, atol=1e-12):
                raise ValueError("Saved predictions and log loss disagree")


def weights(frame: pd.DataFrame, mode: str, *, balance_donors: bool) -> np.ndarray:
    if mode == "clone":
        w = 1 / frame.groupby(["patient", "clone_key"]).clone_key.transform("size")
    elif mode == "cell":
        w = pd.Series(1.0, index=frame.index)
    else:
        raise ValueError(mode)
    if balance_donors:
        w = w / w.groupby(frame.patient).transform("sum")
        w = w * len(frame) / w.sum()
    return w.to_numpy()


def training_support(train: pd.DataFrame, column: str, minimum: int) -> set[str]:
    unique = train.drop_duplicates(["patient", "clone_key", column])
    counts = unique.groupby(column).agg(n=("clone_key", "size"), donors=("patient", "nunique"))
    return set(counts.index[(counts.n >= minimum) & (counts.donors >= 3)])


def load_cached(dataset: Path, reference: Path) -> tuple[pd.DataFrame, dict]:
    manifest_path = dataset / "rfu/run_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    for file, key in (
        ("km5000noMax.Rdata", "km5000_rdata_sha256"),
        ("RFU.R", "rfu_r_sha256"),
        ("trimerMDSfit_small.Rdata", "trimer_rdata_sha256"),
    ):
        if sha256(reference / file) != manifest[key]:
            raise ValueError(f"Cached assignment reference mismatch: {file}")
    if manifest["rfu_threshold"] != 0.6 or manifest["backend_mode"] != "standard":
        raise ValueError("Assignment settings differ from frozen specification")
    tcr_path = dataset / "GSE190905_TCR_data.csv.gz"
    assignment_path = dataset / "rfu/rfu_results_per_row.tsv.gz"
    source = pd.read_csv(tcr_path, index_col=0)
    assigned = pd.read_csv(assignment_path, sep="\t")
    frame = assigned.merge(
        source[["patient", "cluster", "state", "method"]],
        left_on="cell_id",
        right_index=True,
        validate="one_to_one",
    )
    if len(frame) != len(assigned):
        raise ValueError("Missing source cell labels")
    # Upstream R which.max indexes centroids from one. No renumbering is allowed.
    if not (frame.rfu_label == "RFU" + frame.rfu_id.astype(str)).all():
        raise ValueError("RFU label/index numbering mismatch")
    keep = (
        frame.rfu_pass_threshold.eq(True)
        & frame.v_call.fillna("").str.fullmatch(r"TRBV[0-9]+(?:-[0-9]+)?(?:\*[0-9]+)?")
        & frame.j_call.fillna("").str.fullmatch(r"TRBJ[0-9]+-[0-9]+(?:\*[0-9]+)?")
        & frame.cdr3aa.fillna("").str.fullmatch("[ACDEFGHIKLMNPQRSTVWY]+")
    )
    eligible = frame.loc[keep].copy()
    eligible["length"] = eligible.cdr3aa.str.len().astype(str)
    eligible["clone_key"] = eligible[["v_call", "j_call", "cdr3aa"]].agg("|".join, axis=1)
    counts = eligible.patient.value_counts()
    if len(counts) < 4 or counts.min() < 100:
        raise ValueError("Insufficient donor support for frozen comparison")
    info = {
        "source_cells": len(source),
        "assigned_cells": len(frame),
        "eligible_cells": len(eligible),
        "donors": counts.to_dict(),
        "source_clusters": source.cluster.value_counts().to_dict(),
        "input_sha256": {str(p): sha256(p) for p in [tcr_path, assignment_path, manifest_path]},
        "assignment_manifest": manifest,
        "rfuwas_number_compatibility": "not independently checksum-verified",
        "labels": "source expression clusters; source TCR-HVG exclusion undocumented",
    }
    return eligible, info


def save_composition(frame: pd.DataFrame, out: Path) -> None:
    # Descriptive composition uses real V calls; donor clones, not repeated cells.
    composition = (
        frame.drop_duplicates(["patient", "clone_key"])
        .groupby(["rfu_id", "v_call"], as_index=False)
        .agg(clones=("clone_key", "size"), donors=("patient", "nunique"))
    )
    composition["fraction"] = composition.clones / composition.groupby("rfu_id").clones.transform(
        "sum"
    )
    composition.to_csv(out / "observed_rfu_trbv_composition.tsv", sep="\t", index=False)
    unique = frame.drop_duplicates(["patient", "clone_key"])
    summary = unique.groupby("rfu_id", as_index=False).agg(
        donor_clones=("clone_key", "size"),
        donors=("patient", "nunique"),
        distinct_trbv=("v_call", "nunique"),
    )
    dominant = composition.sort_values(
        ["clones", "v_call"], ascending=[False, True]
    ).drop_duplicates("rfu_id")
    summary = summary.merge(
        dominant[["rfu_id", "v_call", "fraction"]], on="rfu_id", validate="one_to_one"
    )
    summary.to_csv(out / "observed_rfu_summary.tsv", sep="\t", index=False)


def run(dataset: Path, reference: Path, out: Path, *, purge_shared: bool = True) -> None:
    out.mkdir(parents=True, exist_ok=True)
    status_path = out / "completion.json"
    status_path.write_text(json.dumps({"status": "running", "purge_shared": purge_shared}) + "\n")
    frame, info = load_cached(dataset, reference)
    save_composition(frame, out)
    rows, coverage, predictions = [], [], []
    for donor in sorted(frame.patient.unique()):
        train, test = split_donor(frame, donor, purge_shared=purge_shared)
        supported_classes = training_support(train, "cluster", 20)
        n_test = len(test)
        train = train.loc[train.cluster.isin(supported_classes)].copy()
        test = test.loc[test.cluster.isin(supported_classes)].copy()
        supported_rfus = training_support(train, "rfu_label", 10)
        for part in (train, test):
            part["rfu_feature"] = part.rfu_label.where(part.rfu_label.isin(supported_rfus), "rare")
        coverage.append(
            {
                "held_out_donor": donor,
                "train_cells": len(train),
                "test_cells": len(test),
                "test_before_class_filter": n_test,
                "classes": sorted(supported_classes),
                "supported_rfus": len(supported_rfus),
                "test_supported_rfu_fraction": float(test.rfu_label.isin(supported_rfus).mean()),
                "train_donors": sorted(train.patient.unique()),
                "train_cell_ids_sha256": hashlib.sha256(
                    "\n".join(sorted(train.cell_id)).encode()
                ).hexdigest(),
                "test_cell_ids_sha256": hashlib.sha256(
                    "\n".join(sorted(test.cell_id)).encode()
                ).hexdigest(),
                "public_clone_purged_cells": int(
                    (
                        (frame.patient != donor)
                        & frame.clone_key.isin(frame.loc[frame.patient == donor, "clone_key"])
                    ).sum()
                )
                if purge_shared
                else 0,
            }
        )
        for mode in ("clone", "cell"):
            train_w = weights(train, mode, balance_donors=True)
            test_w = weights(test, mode, balance_donors=False)
            for model_name in ("TRBV_TRBJ_length", "TRBV_TRBJ_length_RFU"):
                columns = ["v_call", "j_call", "length"]
                if model_name.endswith("_RFU"):
                    columns.append("rfu_feature")
                encoder = OneHotEncoder(handle_unknown="ignore", dtype=np.float64)
                x_train = encoder.fit_transform(train[columns])
                x_test = encoder.transform(test[columns])
                model = LogisticRegression(C=1.0, solver="lbfgs", max_iter=1000, tol=1e-6)
                model.fit(x_train, train.cluster, sample_weight=train_w)
                if model.n_iter_.max() >= model.max_iter:
                    raise ValueError("Model did not converge")
                prob = model.predict_proba(x_test)
                pred = test[["cell_id", "cluster"]].reset_index(drop=True).copy()
                pred["held_out_donor"], pred["weighting"], pred["model"] = donor, mode, model_name
                pred["evaluation_weight"] = test_w
                for i, label in enumerate(model.classes_):
                    pred[f"prob:{label}"] = prob[:, i]
                predictions.append(pred)
                rows.append(
                    {
                        "held_out_donor": donor,
                        "weighting": mode,
                        "model": model_name,
                        "test_cells": len(test),
                        "test_donor_clones": test.clone_key.nunique(),
                        "log_loss": log_loss(
                            test.cluster, prob, labels=model.classes_, sample_weight=test_w
                        ),
                        "balanced_accuracy": balanced_accuracy_score(
                            test.cluster, model.classes_[prob.argmax(axis=1)], sample_weight=test_w
                        ),
                    }
                )
        print(f"Finished held-out donor {donor}", flush=True)
    result = pd.DataFrame(rows)
    predicted = pd.concat(predictions, ignore_index=True)
    audit_predictions(predicted, result)
    predicted.to_csv(out / "heldout_predictions.tsv.gz", sep="\t", index=False)
    paired_deltas(result).to_csv(out / "paired_donor_deltas.tsv", sep="\t", index=False)
    result.to_csv(out / "donor_held_out_metrics.tsv", sep="\t", index=False)
    info["folds"] = coverage
    info["settings"] = {
        "C": 1.0,
        "threshold": 0.6,
        "minimum_class_clones": 20,
        "minimum_rfu_clones": 10,
        "minimum_training_donors": 3,
        "public_receptor_purge": purge_shared,
        "estimand": "unseen_receptor_donor_generalization"
        if purge_shared
        else "donor_generalization",
        "hyperparameter_search": False,
        "clone_definition": "donor + primary TRBV/TRBJ/CDR3aa, across timepoints",
    }
    info["python"] = sys.executable
    info["git_sha"] = subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    info["script_sha256"] = sha256(Path(__file__))
    (out / "provenance.json").write_text(json.dumps(info, indent=2) + "\n")
    outputs = [
        "heldout_predictions.tsv.gz",
        "paired_donor_deltas.tsv",
        "donor_held_out_metrics.tsv",
        "provenance.json",
        "observed_rfu_trbv_composition.tsv",
        "observed_rfu_summary.tsv",
    ]
    status_path.write_text(
        json.dumps(
            {
                "status": "complete",
                "predictions_audited": True,
                "script_sha256": sha256(Path(__file__)),
                "outputs": {name: sha256(out / name) for name in outputs},
            },
            indent=2,
        )
        + "\n"
    )
    print(
        result.groupby(["weighting", "model"])[["log_loss", "balanced_accuracy"]].mean().to_string()
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument(
        "--retain-shared-receptors",
        action="store_true",
        help="Donor-generalization sensitivity; original default excludes shared receptors",
    )
    args = parser.parse_args()
    run(args.dataset, args.reference, args.out, purge_shared=not args.retain_shared_receptors)
