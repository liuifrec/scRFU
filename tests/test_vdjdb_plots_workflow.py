from __future__ import annotations

import json

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from matplotlib.axes import Axes

from examples.vdjdb_antigen_evidence import build_parser, run_analysis
from scrfu.pl import (
    antigen_permutation_distribution,
    rfu_antigen_bubble,
    rfu_antigen_heatmap,
)
from scrfu.pl import rfu_antigen_coherence as plot_rfu_antigen_coherence
from scrfu.tl import AntigenPermutationResult

matplotlib.use("Agg")


def _abundance() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "rfu_label": ["R1", "R1", "R2"],
            "epitope": ["A", "B", "A"],
            "antigen_abundance": [3.0, 1.0, 2.0],
            "within_rfu_proportion": [0.75, 0.25, 1.0],
        }
    )


def test_antigen_plots_return_axes_without_show(monkeypatch) -> None:
    monkeypatch.setattr(plt, "show", lambda: pytest.fail("plt.show() must not be called"))
    coherence = pd.DataFrame(
        {
            "rfu_label": ["R1", "R2"],
            "vdjdb_matched_sequences": [3, 2],
            "antigen_purity": [0.75, 1.0],
            "total_rfu_sequences": [4, 2],
            "eligible_for_coherence": [True, True],
        }
    )
    permutation = AntigenPermutationResult(
        observed=0.8,
        permutation_values=np.array([0.1, 0.3, 0.5]),
        null_mean=0.3,
        null_std=0.2,
        empirical_upper_tail_probability=0.25,
        z_score=2.5,
        parameters={"metric": "same_antigen_pair_fraction"},
    )
    axes = [
        rfu_antigen_heatmap(_abundance(), top_n_rfus=1),
        plot_rfu_antigen_coherence(coherence, color="eligible_for_coherence"),
        antigen_permutation_distribution(permutation, bins=3),
        rfu_antigen_bubble(_abundance(), top_n_antigens=1),
    ]
    assert all(isinstance(ax, Axes) for ax in axes)
    assert axes[0].get_xlabel() == "epitope"
    assert "Permutation" in axes[2].get_ylabel()


@pytest.mark.parametrize(
    "function,args",
    [
        (rfu_antigen_heatmap, (pd.DataFrame(),)),
        (rfu_antigen_bubble, (pd.DataFrame(),)),
        (
            plot_rfu_antigen_coherence,
            (
                pd.DataFrame(
                    columns=[
                        "vdjdb_matched_sequences",
                        "antigen_purity",
                        "total_rfu_sequences",
                    ]
                ),
            ),
        ),
    ],
)
def test_antigen_plots_reject_empty(function, args) -> None:
    with pytest.raises(ValueError, match="empty|no finite"):
        function(*args)


def _workflow_sequences() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "input_row_id": ["q1", "q2", "q3", "q4"],
            "unique_sequence_id": ["s1", "s2", "s3", "s4"],
            "cdr3aa": ["CASSAAA", "CASSAAD", "CASSCCC", "CASSCCD"],
            "v_call": ["TRBV1", "TRBV1", "TRBV2", "TRBV2"],
            "chain": ["TRB"] * 4,
            "rfu_label": ["R1", "R1", "R2", "R2"],
            "pass_thr": [True, True, True, True],
        }
    )


def test_vdjdb_workflow_help() -> None:
    with pytest.raises(SystemExit) as error:
        build_parser().parse_args(["--help"])
    assert error.value.code == 0


def test_vdjdb_workflow_offline_end_to_end(tmp_path) -> None:
    sequence_path = tmp_path / "rfu_sequences.tsv.gz"
    row_path = tmp_path / "rfu_rows.tsv.gz"
    reference_path = tmp_path / "vdjdb.tsv.gz"
    metadata_path = tmp_path / "metadata.tsv.gz"
    outdir = tmp_path / "out"
    sequences = _workflow_sequences()
    sequences.to_csv(sequence_path, sep="\t", index=False)
    rows = pd.concat(
        [
            sequences.assign(cell_id=["c1", "c2", "c3", "c4"]),
            sequences.iloc[[0]].assign(input_row_id="q5", cell_id="c5"),
        ],
        ignore_index=True,
    )
    rows.to_csv(row_path, sep="\t", index=False)
    pd.DataFrame(
        {
            "cdr3aa": sequences["cdr3aa"],
            "v_call": sequences["v_call"],
            "chain": ["TRB"] * 4,
            "epitope": ["A", "A", "A", "B"],
            "score": [3, 2, 2, 1],
        }
    ).to_csv(reference_path, sep="\t", index=False)
    pd.DataFrame(
        {
            "cell_id": ["c1", "c2", "c3", "c4", "c5"],
            "sample": ["p1", "p1", "p2", "p2", "p2"],
            "donor": ["d1", "d1", "d2", "d2", "d2"],
            "phenotype": ["T1", "T2", "T1", "T2", "T1"],
        }
    ).to_csv(metadata_path, sep="\t", index=False)
    manifest = run_analysis(
        rfu_sequences=sequence_path,
        rfu_rows=row_path,
        vdjdb=reference_path,
        vdjdb_release="synthetic-release",
        outdir=outdir,
        metadata=metadata_path,
        sample_key="sample",
        donor_key="donor",
        phenotype_key="phenotype",
        n_permutations=10,
        random_state=4,
        save_permutation_values=True,
    )
    expected = {
        "vdjdb_reference_qc.json",
        "vdjdb_matches_long.tsv.gz",
        "sequence_antigen_summary.tsv.gz",
        "row_antigen_summary.tsv.gz",
        "rfu_antigen_coherence.tsv.gz",
        "rfu_antigen_abundance.tsv.gz",
        "global_antigen_coherence.json",
        "antigen_grouping_comparison.tsv.gz",
        "permutation_summary.json",
        "permutation_values.tsv.gz",
        "rfu_antigen_heatmap.png",
        "rfu_antigen_coherence.png",
        "antigen_permutation_distribution.png",
        "rfu_antigen_bubble.png",
        "sequence_sample_prevalence.tsv.gz",
        "rfu_antigen_recurrence.tsv.gz",
        "rfu_antigen_phenotype_abundance.tsv.gz",
        "antigen_evidence_tiers_by_phenotype.tsv.gz",
        "antigen_ambiguity_by_phenotype.tsv.gz",
        "run_manifest.json",
    }
    assert expected == {path.name for path in outdir.iterdir()}
    assert manifest["dimensions"]["matched_unique_sequences"] == 4
    assert manifest["join_qc"]["unmatched_result_cells"] == 0
    saved = json.loads((outdir / "run_manifest.json").read_text())
    assert saved["parameters"]["vdjdb_release"] == "synthetic-release"
    assert saved["permutation"]["status"] == "ok"
