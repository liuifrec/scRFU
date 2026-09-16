"""Synthetic public-resource selection and paired-barcode coverage checks."""

import gzip

import pandas as pd
import pytest

from manuscript.scripts.gse280982_prepare import (
    obtain,
    pair_receptors,
    select_files,
    treatment_label,
    validate_gzip,
)


def resources():
    prefix = "https://ftp.ncbi.nlm.nih.gov/geo/samples/synthetic/"
    inventory = pd.DataFrame(
        {
            "donor": [1, 1, 1, 2],
            "visit": [1] * 4,
            "compartment": ["tumor"] * 4,
            "assay": ["TCR", "GEX", "GEX", "TCR"],
            "longitudinal_hypr_sample": [True, True, True, False],
            "url": [
                prefix + x
                for x in [
                    "filtered_contig_annotations.csv.gz",
                    "barcodes.tsv.gz",
                    "matrix.mtx.gz",
                    "other_filtered_contig_annotations.csv.gz",
                ]
            ],
        }
    )
    pairs = pd.DataFrame(
        {
            "donor": [1],
            "visit": [1],
            "compartment": ["tumor"],
            "paired_processed_GEX_TCR_available": [True],
        }
    )
    return inventory, pairs


def test_download_selection_excludes_matrices_and_other_cohorts():
    inventory, pairs = resources()
    selected = select_files(inventory, pairs)
    assert len(selected) == 2 and set(selected.assay) == {"TCR", "GEX"}
    with pytest.raises(ValueError, match="Multiple"):
        select_files(pd.concat([inventory, inventory.iloc[[0]]]), pairs)
    inventory.loc[0, "url"] = "https://unapproved.example/filtered_contig_annotations.csv.gz"
    with pytest.raises(ValueError, match="HTTPS NCBI"):
        select_files(inventory, pairs)


def test_source_treatment_defines_time_not_filename_alone():
    assert (
        treatment_label("tissue: cancer; treatment: Last Day of Radition", 2)
        == "Last Day of Radiation"
    )
    assert treatment_label("treatment: 6 Weeks Post-Radiation", 3) == "6 Weeks Post-Radiation"
    with pytest.raises(ValueError, match="Visit number"):
        treatment_label("treatment: Pre-Treatment", 2)


def test_primary_chain_selection_then_gex_matching_preserves_missing_cells():
    contigs = pd.DataFrame(
        {
            "barcode": ["a", "a", "b", "c"],
            "chain": ["TRB"] * 4,
            "cdr3": ["CASSF", "CATSF", "CASRF", "CASKF"],
            "v_gene": ["TRBV1"] * 4,
            "j_gene": ["TRBJ1-1"] * 4,
            "productive": [True, True, True, False],
            "is_cell": [True] * 4,
            "high_confidence": [True] * 4,
            "umis": [2, 8, 5, 9],
            "reads": [20, 80, 50, 90],
        }
    )
    selected, qc = pair_receptors(contigs, pd.Series(["a", "d"]), "synthetic_visit1")
    assert selected.cdr3aa.tolist() == ["CATSF"]
    assert selected.cell_id.tolist() == ["synthetic_visit1:a"]
    assert qc["primary_TRB_cells"] == 2
    assert qc["primary_TRB_cells_in_GEX"] == 1 and qc["primary_TRB_cells_not_in_GEX"] == 1
    assert qc["cells_with_multiple_productive_TRB"] == 1
    assert qc["gex_barcodes"] == 2
    with pytest.raises(ValueError, match="unique"):
        pair_receptors(contigs, pd.Series(["a", "a"]), "synthetic")


def test_absent_resource_stays_unknown_and_gzip_integrity_is_checked(tmp_path):
    result = obtain("https://ftp.ncbi.nlm.nih.gov/synthetic.csv.gz", tmp_path, False)
    assert result["status"] == "missing_not_downloaded"
    assert "bytes" not in result  # unknown is not a zero-byte observed repertoire
    path = tmp_path / "synthetic.csv.gz"
    with gzip.open(path, "wt") as stream:
        stream.write("synthetic\n1\n")
    validate_gzip(path)
    data = bytearray(path.read_bytes())
    data[-5] ^= 0xFF
    path.write_bytes(data)
    with pytest.raises((OSError, EOFError)):
        validate_gzip(path)
