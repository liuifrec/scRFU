from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import anndata as ad
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "examples" / "gse190905_prepare.py"


def _subprocess_output(result: subprocess.CompletedProcess[str]) -> str:
    return f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"


def _write_csv_gz(path: Path, df: pd.DataFrame) -> None:
    df.to_csv(path, index=False, compression="gzip")


def test_gse190905_prepare_automatic_column_inference(tmp_path: Path) -> None:
    tcr_path = tmp_path / "tcr.csv.gz"
    metadata_path = tmp_path / "meta.csv.gz"
    out_path = tmp_path / "gse190905_scrfu_input.h5ad"
    report_path = tmp_path / "report.json"

    _write_csv_gz(
        tcr_path,
        pd.DataFrame(
            {
                "barcode": ["c1", "c2", "x1"],
                "locus": ["TRB", "TRB", "TRA"],
                "junction_aa": ["CASSA", "CASST", "CAVR"],
                "v_gene": ["TRBV1", "TRBV2", "TRAV1"],
                "productive": [True, False, True],
            }
        ),
    )
    _write_csv_gz(
        metadata_path,
        pd.DataFrame({"cell_id": ["c1", "c2", "c3"], "sample": ["s1", "s1", "s2"]}),
    )

    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--tcr",
            str(tcr_path),
            "--metadata",
            str(metadata_path),
            "--out",
            str(out_path),
            "--report",
            str(report_path),
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, _subprocess_output(result)
    assert out_path.exists()
    assert report_path.exists()

    report = json.loads(report_path.read_text(encoding="utf-8"))
    assert report["n_metadata_rows"] == 3
    assert report["n_tcr_rows"] == 3
    assert report["n_cells_with_tcr"] == 3
    assert report["n_overlap_cells"] == 2
    assert report["barcode_overlap_rate"] == 2 / 3
    assert report["inferred_columns"]["cdr3_col"] == "junction_aa"
    assert report["n_productive_trb_rows"] == 1

    adata = ad.read_h5ad(out_path)
    assert "airr" in adata.obsm
    assert list(adata.obs_names) == ["c1", "c2", "c3"]
    assert list(adata.obsm["airr"]["cell_id"]) == ["c1", "c2", "c3"]


def test_gse190905_prepare_explicit_column_mapping(tmp_path: Path) -> None:
    tcr_path = tmp_path / "tcr.csv.gz"
    metadata_path = tmp_path / "meta.csv.gz"
    out_path = tmp_path / "prepared.h5ad"
    report_path = tmp_path / "report.json"

    _write_csv_gz(
        tcr_path,
        pd.DataFrame(
            {
                "CellBarcode": ["c1"],
                "ChainName": ["TRB"],
                "AASeq": ["CASSA"],
                "VGene": ["TRBV1"],
                "IsProductive": ["yes"],
            }
        ),
    )
    _write_csv_gz(metadata_path, pd.DataFrame({"BarcodeID": ["c1"], "sample": ["s1"]}))

    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--tcr",
            str(tcr_path),
            "--metadata",
            str(metadata_path),
            "--out",
            str(out_path),
            "--report",
            str(report_path),
            "--cell-col",
            "CellBarcode",
            "--chain-col",
            "ChainName",
            "--cdr3-col",
            "AASeq",
            "--v-col",
            "VGene",
            "--productive-col",
            "IsProductive",
            "--metadata-cell-col",
            "BarcodeID",
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, _subprocess_output(result)
    report = json.loads(report_path.read_text(encoding="utf-8"))
    assert report["inferred_columns"]["tcr_cell_col"] == "CellBarcode"
    assert report["inferred_columns"]["metadata_cell_col"] == "BarcodeID"
    assert report["n_productive_trb_rows"] == 1


def test_gse190905_prepare_missing_cdr3_column_fails_clearly(tmp_path: Path) -> None:
    tcr_path = tmp_path / "tcr.csv.gz"
    metadata_path = tmp_path / "meta.csv.gz"

    _write_csv_gz(
        tcr_path,
        pd.DataFrame({"cell_id": ["c1"], "chain": ["TRB"], "v_call": ["TRBV1"]}),
    )
    _write_csv_gz(metadata_path, pd.DataFrame({"cell_id": ["c1"]}))

    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--tcr",
            str(tcr_path),
            "--metadata",
            str(metadata_path),
            "--out",
            str(tmp_path / "out.h5ad"),
            "--report",
            str(tmp_path / "report.json"),
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )

    assert result.returncode != 0, _subprocess_output(result)
    assert "could not infer required TCR column for cdr3" in result.stderr, _subprocess_output(
        result
    )
