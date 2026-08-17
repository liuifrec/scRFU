from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest
from anndata import AnnData

from scrfu.cli import main as cli_main
from scrfu.io import write_receptor_cache
from scrfu.pp import canonicalize_receptor_table
from scrfu.rfu import RFURunResult
from scrfu.tl import call_rfu, call_rfu_table, rfu_metrics

ROOT = Path(__file__).resolve().parents[1]
WORKFLOW = ROOT / "examples" / "receptor_table_workflow.py"


def _rfu_dir(path: Path) -> Path:
    path.mkdir()
    (path / "RFU.R").write_text("AssignRFUs <- function(ff) {}\n", encoding="utf-8")
    (path / "trimerMDSfit_small.Rdata").write_bytes(b"trimer")
    (path / "km5000noMax.Rdata").write_bytes(b"centers")
    return path


def _wrapper(path: Path) -> Path:
    path.write_text("#!/usr/bin/env Rscript\n", encoding="utf-8")
    return path


def _canonical(adapter: str = "synthetic") -> pd.DataFrame:
    return canonicalize_receptor_table(
        pd.DataFrame(
            {
                "cell_id": ["c1", "c2", "c3", "c4"],
                "chain": ["TRB", "TRB", "TRB", "TRA"],
                "cdr3aa": ["CASSA", "CASSA", "ASSBAD", "CAVR"],
                "v_call": ["TRBV1", "TRBV2", "TRBV3", "TRAV1"],
                "productive": True,
                "source_adapter": adapter,
                "source_row_id": ["a", "b", "c", "d"],
            }
        )
    )


def _fake_run(self: object, features: pd.DataFrame, **kwargs: object) -> RFURunResult:
    del self, kwargs
    rows = features.copy().reset_index(drop=True)
    eligible = rows["cdr3aa"].str.startswith("C")
    rows["eligibility_status"] = "ineligible_cdr3_not_starting_c"
    rows.loc[eligible, "eligibility_status"] = "eligible"
    rows["unique_sequence_id"] = pd.Series(pd.NA, index=rows.index, dtype="string")
    rows.loc[eligible, "unique_sequence_id"] = "sequence_00000000"
    rows["rfu_id"] = pd.NA
    rows.loc[eligible, "rfu_id"] = 7
    rows["rfu_label"] = pd.NA
    rows.loc[eligible, "rfu_label"] = "RFU7"
    rows["rfu_score"] = pd.NA
    rows.loc[eligible, "rfu_score"] = 0.9
    rows["pass_thr"] = pd.Series(pd.NA, index=rows.index, dtype="boolean")
    rows.loc[eligible, "pass_thr"] = True
    rows["rfu_status"] = rows["eligibility_status"]
    rows.loc[eligible, "rfu_status"] = "assigned_threshold_pass"
    return RFURunResult(
        rows,
        "",
        "",
        0,
        metadata={
            "original_row_count": len(rows),
            "eligible_row_count": int(eligible.sum()),
            "unique_query_count": 1,
            "reconstructed_output_row_count": len(rows),
        },
    )


def test_call_rfu_table_preserves_canonical_provenance_and_reconstructs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr("scrfu.backends.rfu_repo.RFURepoBackend.run", _fake_run)
    result = call_rfu_table(
        _canonical(),
        rfu_dir=_rfu_dir(tmp_path / "rfu"),
        wrapper_r_path=_wrapper(tmp_path / "wrapper.R"),
    )

    assert result.per_row["input_row_id"].tolist() == [
        "row_00000000",
        "row_00000001",
        "row_00000002",
    ]
    assert result.per_row["source_row_id"].tolist() == ["a", "b", "c"]
    assert result.per_sequence["multiplicity"].tolist() == [2]
    assert result.mapping["unique_sequence_id"].tolist()[:2] == [
        "sequence_00000000",
        "sequence_00000000",
    ]
    assert result.per_row["rfu_id_nearest"].equals(result.per_row["rfu_id"])
    assert result.per_row["rfu_pass_threshold"].equals(result.per_row["pass_thr"].astype("boolean"))
    assert result.per_row["reference_coverage_status"].tolist() == [
        "covered",
        "covered",
        "ineligible_sequence",
    ]
    assert result.provenance["table_level_api"] is True


def test_equivalent_adapter_content_has_identical_rfu_and_metrics(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr("scrfu.backends.rfu_repo.RFURepoBackend.run", _fake_run)
    rfu_dir = _rfu_dir(tmp_path / "rfu")
    wrapper = _wrapper(tmp_path / "wrapper.R")
    left = call_rfu_table(_canonical("wells_tcr_ir"), rfu_dir=rfu_dir, wrapper_r_path=wrapper)
    right = call_rfu_table(_canonical("scirpy_airr"), rfu_dir=rfu_dir, wrapper_r_path=wrapper)
    columns = ["cell_id", "cdr3aa", "rfu_id", "rfu_label", "rfu_score", "pass_thr"]
    pd.testing.assert_frame_equal(left.per_row[columns], right.per_row[columns])

    for result in (left, right):
        result.per_row["phenotype"] = "p"
    left_metrics = rfu_metrics(left.per_row, groupby="phenotype", weighting="cell", chain="TRB")
    right_metrics = rfu_metrics(right.per_row, groupby="phenotype", weighting="cell", chain="TRB")
    pd.testing.assert_frame_equal(left_metrics, right_metrics)


def test_existing_anndata_call_delegates_with_equivalent_results(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr("scrfu.backends.rfu_repo.RFURepoBackend.run", _fake_run)
    rfu_dir = _rfu_dir(tmp_path / "rfu")
    wrapper = _wrapper(tmp_path / "wrapper.R")
    airr = pd.DataFrame(
        {
            "cell_id": ["c1", "c2", "c3"],
            "chain": ["TRB", "TRB", "TRB"],
            "cdr3aa": ["CASSA", "CASSA", "ASSBAD"],
            "v_call": ["TRBV1", "TRBV2", "TRBV3"],
            "productive": True,
        },
        index=["c1", "c2", "c3"],
    )
    adata = AnnData(obs=pd.DataFrame(index=["c1", "c2", "c3"]))
    adata.obsm["airr"] = airr
    ann_result = call_rfu(adata, rfu_dir=rfu_dir, wrapper_r_path=wrapper)
    table_result = call_rfu_table(
        _canonical().iloc[:3], rfu_dir=rfu_dir, wrapper_r_path=wrapper
    ).per_row

    columns = ["cell_id", "cdr3aa", "rfu_id", "rfu_label", "rfu_score", "pass_thr"]
    pd.testing.assert_frame_equal(
        ann_result[columns].reset_index(drop=True), table_result[columns].reset_index(drop=True)
    )
    assert adata.uns["scrfu"]["rfu"]["table_level_api"] is True


def test_metrics_assignment_policy_is_explicit() -> None:
    frame = pd.DataFrame(
        {
            "cell_id": ["c1", "c2"],
            "chain": ["TRB", "TRB"],
            "cdr3aa": ["CASSA", "CASST"],
            "rfu_label": ["RFU1", "RFU1"],
            "pass_thr": [True, False],
            "phenotype": ["p", "p"],
        }
    )
    nearest = rfu_metrics(
        frame, groupby="phenotype", weighting="cell", assignment_policy="nearest", chain="TRB"
    )
    passed = rfu_metrics(
        frame,
        groupby="phenotype",
        weighting="cell",
        assignment_policy="threshold_pass",
        chain="TRB",
    )
    assert nearest.loc[0, "rfu_cell_count"] == 2
    assert nearest.loc[0, "rfu_threshold_pass_rate"] == 0.5
    assert passed.loc[0, "rfu_cell_count"] == 1
    assert passed.loc[0, "rfu_threshold_pass_rate"] == 1.0


def test_prepare_receptors_cli_lists_adapters_and_writes_dataframe_cache(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    cli_main(["prepare-receptors", "--list-adapters"])
    assert "wells_tcr_ir" in capsys.readouterr().out

    table = pd.DataFrame(
        {
            "cell_id": ["c1"],
            "chain": ["TRB"],
            "cdr3aa": ["CASSA"],
            "v_call": ["TRBV1"],
            "productive": [True],
        }
    )
    source = tmp_path / "airr.tsv"
    table.to_csv(source, sep="\t", index=False)
    cache = tmp_path / "cache"
    cli_main(
        [
            "prepare-receptors",
            "--input",
            str(source),
            "--adapter",
            "generic_airr_dataframe",
            "--outdir",
            str(cache),
        ]
    )
    assert (cache / "receptors.tsv.gz").is_file()
    assert (cache / "preparation_manifest.json").is_file()


def test_generic_workflow_help_and_offline_canonical_table(tmp_path: Path) -> None:
    help_result = subprocess.run(
        [sys.executable, str(WORKFLOW), "--help"], capture_output=True, text=True
    )
    assert help_result.returncode == 0
    assert "--cache" in help_result.stdout

    receptors = tmp_path / "receptors.tsv"
    _canonical().to_csv(receptors, sep="\t", index=False)
    outdir = tmp_path / "out"
    result = subprocess.run(
        [
            sys.executable,
            str(WORKFLOW),
            "--receptors",
            str(receptors),
            "--outdir",
            str(outdir),
            "--skip-rfu",
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    assert {
        "receptors.tsv.gz",
        "unique_sequence_map.tsv.gz",
        "rfu_results_per_sequence.tsv.gz",
        "rfu_results_per_row.tsv.gz",
        "run_manifest.json",
    }.issubset(path.name for path in outdir.iterdir())


def test_generic_workflow_offline_cache_input(tmp_path: Path) -> None:
    cache = tmp_path / "cache"
    write_receptor_cache(
        cache,
        _canonical(),
        pd.DataFrame(index=["c1", "c2", "c3", "c4"]),
        source_adapter="synthetic",
        source_adapter_version="1",
        source_format="table",
    )
    outdir = tmp_path / "out"
    result = subprocess.run(
        [
            sys.executable,
            str(WORKFLOW),
            "--cache",
            str(cache),
            "--outdir",
            str(outdir),
            "--skip-rfu",
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    manifest = (outdir / "run_manifest.json").read_text(encoding="utf-8")
    assert '"input_kind": "receptor_cache"' in manifest
    assert '"cache_schema_version": 2' in manifest


def test_generic_workflow_offline_scirpy_h5ad_input(tmp_path: Path) -> None:
    adata = AnnData(obs=pd.DataFrame({"phenotype": ["p", "p"]}, index=["c1", "c2"]))
    adata.obsm["airr"] = pd.DataFrame(
        {
            "cell_id": ["c1", "c2"],
            "chain": ["TRB", "TRB"],
            "cdr3aa": ["CASSA", "CASST"],
            "v_call": ["TRBV1", "TRBV2"],
            "productive": [True, True],
        },
        index=adata.obs_names,
    )
    source = tmp_path / "airr.h5ad"
    adata.write_h5ad(source)
    outdir = tmp_path / "out"
    result = subprocess.run(
        [
            sys.executable,
            str(WORKFLOW),
            "--input",
            str(source),
            "--adapter",
            "scirpy_airr",
            "--metadata-column",
            "phenotype",
            "--outdir",
            str(outdir),
            "--skip-rfu",
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    assert (outdir / "adapter_qc.json").is_file()
    assert len(pd.read_csv(outdir / "rfu_results_per_row.tsv.gz", sep="\t")) == 2


def test_original_rfu_table_backend_rejects_non_trb_before_execution() -> None:
    rows = pd.DataFrame(
        {
            "input_row_id": ["r1"],
            "cell_id": ["c1"],
            "chain": ["IGH"],
            "cdr3aa": ["CARF"],
            "v_call": ["IGHV1"],
            "source_adapter": ["synthetic"],
            "source_row_id": ["source-1"],
        }
    )
    with pytest.raises(ValueError, match="TRB-specific"):
        call_rfu_table(rows, chain="IGH")
