from __future__ import annotations

import pandas as pd
import pytest

from scrfu.adapters import adapt_cellranger_vdj, get_receptor_adapter
from scrfu.cli import main


def _contigs() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "barcode": ["AA-1", "AA-1", "BB-1", "CC-1", "DD-1"],
            "is_cell": [True, True, False, True, True],
            "contig_id": ["a_low", "a_high", "b", "c", "d"],
            "high_confidence": [True, True, True, False, True],
            "chain": ["TRB", "TRB", "TRB", "TRA", "IGH"],
            "v_gene": ["TRBV1", "TRBV2", "TRBV3", "TRAV1", "IGHV1"],
            "d_gene": ["TRBD1", "TRBD1", "TRBD2", None, "IGHD1"],
            "j_gene": ["TRBJ1", "TRBJ2", "TRBJ2", "TRAJ1", "IGHJ1"],
            "c_gene": ["TRBC1", "TRBC2", "TRBC1", "TRAC", "IGHM"],
            "productive": [True, True, True, True, False],
            "cdr3": ["CASSL", "CASSH", "CASSB", "CAVA", "CARW"],
            "cdr3_nt": ["aaa", "bbb", "ccc", "ddd", "eee"],
            "reads": [100, 50, 20, 10, 5],
            "umis": [1, 9, 2, 3, 4],
            "raw_clonotype_id": ["cl1", "cl1", "cl2", "cl3", "cl4"],
            "raw_consensus_id": ["con1", "con2", "con3", "con4", "con5"],
        }
    )


def test_cellranger_primary_priority_and_canonical_fields() -> None:
    result = adapt_cellranger_vdj(_contigs(), chain="TRB")

    assert result.receptors["cell_id"].tolist() == ["AA-1"]
    assert result.receptors["cdr3aa"].tolist() == ["CASSH"]
    assert result.receptors["contig_id"].tolist() == ["a_high"]
    assert result.receptors["junction"].tolist() == ["bbb"]
    assert result.receptors["c_call"].tolist() == ["TRBC2"]
    assert result.receptors["clonotype_id"].tolist() == ["cl1"]
    assert result.receptors["source_row_id"].tolist() == ["1"]


def test_cellranger_filters_all_chains_and_mixed_receptors() -> None:
    all_rows = adapt_cellranger_vdj(
        _contigs(),
        chain=None,
        productive_only=False,
        primary_chain=False,
        filter_is_cell=False,
        filter_high_confidence=False,
    )
    assert all_rows.receptors["chain"].tolist() == ["TRB", "TRB", "TRB", "TRA", "IGH"]

    filtered = adapt_cellranger_vdj(
        _contigs(), chain=None, primary_chain=False, productive_only=True
    )
    assert filtered.receptors["contig_id"].tolist() == ["a_low", "a_high"]


def test_cellranger_missing_optional_fields_and_required_errors() -> None:
    minimal = pd.DataFrame({"barcode": ["x"], "chain": ["beta"], "cdr3": ["CASSF"]})
    result = adapt_cellranger_vdj(minimal)
    assert result.receptors.loc[0, "chain"] == "TRB"
    assert pd.isna(result.receptors.loc[0, "v_call"])

    for missing in ("barcode", "chain", "cdr3"):
        with pytest.raises(ValueError, match="missing required columns"):
            adapt_cellranger_vdj(minimal.drop(columns=missing))


def test_cellranger_csv_dataframe_parity_and_aliases(tmp_path) -> None:
    source = _contigs()
    path = tmp_path / "filtered_contig_annotations.csv"
    source.to_csv(path, index=False)
    frame_result = adapt_cellranger_vdj(source, chain=None, primary_chain=False)
    csv_result = adapt_cellranger_vdj(path, chain=None, primary_chain=False)
    pd.testing.assert_frame_equal(frame_result.receptors, csv_result.receptors)
    assert get_receptor_adapter("cellranger") is adapt_cellranger_vdj
    assert get_receptor_adapter("tenx_vdj") is adapt_cellranger_vdj


def test_cellranger_cli_prepare_receptors(tmp_path, capsys) -> None:
    source = tmp_path / "all_contig_annotations.csv"
    _contigs().to_csv(source, index=False)
    out = tmp_path / "cache"
    main(
        [
            "prepare-receptors",
            "--input",
            str(source),
            "--adapter",
            "cellranger_vdj",
            "--outdir",
            str(out),
            "--chain",
            "TRB",
        ]
    )
    assert (out / "receptors.tsv.gz").is_file()
    assert '"adapter": "cellranger_vdj"' in capsys.readouterr().out
