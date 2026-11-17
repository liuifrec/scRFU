from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

from examples.rfuwas_regulatory_prepare import prepare_rfuwas_data1, read_data1, run
from scrfu.tl import regulatory_triangulation


@pytest.fixture
def data1():
    # Synthetic coordinates/effects, using the verified publication column structure.
    return pd.DataFrame(
        {
            "SNP": ["7_12345_a_g", "6_23456_AT_A"],
            "RFU": [1, 1415],
            "beta": [-0.2, 0.3],
            "t.stat": [-4.0, 5.0],
            "p.value": [1e-8, 1e-9],
        }
    )


def test_canonical_identity_and_unresolved_direction(data1):
    before = data1.copy(deep=True)
    result = prepare_rfuwas_data1(data1, release="synthetic-v1").set_index("SNP")
    assert result.loc["7_12345_a_g", "variant_key"] == "7:12345:A:G"
    assert result.loc["7_12345_a_g", "variant_id"] == "7:12345:A:G"
    assert result.loc["7_12345_a_g", "rfu_label"] == "1"
    assert result.loc["6_23456_AT_A", "rfu_label"] == "1415"
    assert result.loc["7_12345_a_g", "beta"] == -0.2
    assert result.loc["7_12345_a_g", "t.stat"] == -4.0
    assert result.genome_build.eq("GRCh38").all()
    assert result.effect_allele.isna().all() and result.se.isna().all()
    assert result.effect_allele_status.eq("unresolved_in_published_data1").all()
    eqtl = pd.DataFrame(
        {
            "variant_id": ["7:12345:A:G"],
            "genome_build": ["GRCh38"],
            "gene": ["GENE"],
            "source": ["synthetic"],
            "effect_allele": ["G"],
            "beta": [1.0],
        }
    )
    hits = regulatory_triangulation(result.reset_index(), eqtl=eqtl).matched_evidence
    assert len(hits) == 1 and hits.direction_concordant.isna().all()
    assert_frame_equal(data1, before)


@pytest.mark.parametrize("snp", ["rs123", "7_12345_A", "7_12345_A_G,T", "7_0_A_G"])
def test_invalid_variants(data1, snp):
    data1.loc[0, "SNP"] = snp
    with pytest.raises(ValueError):
        prepare_rfuwas_data1(data1, release="v1")


@pytest.mark.parametrize("label", [0, -1, 1.5, None, "RFU1"])
def test_invalid_labels(data1, label):
    data1["RFU"] = data1.RFU.astype(object)
    data1.loc[0, "RFU"] = label
    with pytest.raises(ValueError, match="one-based"):
        prepare_rfuwas_data1(data1, release="v1")


def test_wrong_artifact_and_preexisting_orientation(data1):
    for raw in (
        pd.DataFrame({"REF": ["A"], "ALT": ["G"], "weight": [0.1]}),
        pd.DataFrame({"RFU": [1], "effect": [0.1], "phecode": [1]}),
    ):
        with pytest.raises(ValueError, match="Data 1"):
            prepare_rfuwas_data1(raw, release="v1")
    with pytest.raises(ValueError, match="original Data 1"):
        prepare_rfuwas_data1(data1.assign(effect_allele="G"), release="v1")


@pytest.mark.parametrize(("format", "sep"), [("tsv", "\t"), ("csv", ",")])
def test_file_export_provenance_and_no_overwrite(data1, tmp_path, format, sep):
    path = tmp_path / f"data1.{format}"
    data1.to_csv(path, sep=sep, index=False)
    for directory in ("a", "b"):
        run(path, tmp_path / directory, format=format, release="synthetic", nrows=1)
    for name in ("rfu_qtl.tsv", "provenance.json"):
        assert (tmp_path / "a" / name).read_bytes() == (tmp_path / "b" / name).read_bytes()
    provenance = json.loads((tmp_path / "a" / "provenance.json").read_text())
    assert provenance["rows_written"] == 1 and provenance["rfu_label_offset"] == 0
    assert len(provenance["input_sha256"]) == 64
    with pytest.raises(FileExistsError):
        run(path, tmp_path / "a", format=format, release="v1")
    with pytest.raises(ValueError, match="nrows"):
        read_data1(path, format=format, nrows=0)


def test_xlsx_selects_data1_not_other_sheets(data1, tmp_path):
    pytest.importorskip("openpyxl")
    path = tmp_path / "supplement.xlsx"
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        data1.assign(beta=99).to_excel(writer, sheet_name="Data 4", index=False, startrow=1)
        data1.to_excel(writer, sheet_name="Data 1", index=False, startrow=1)
    result = prepare_rfuwas_data1(read_data1(path, format="xlsx", nrows=1), release="synthetic")
    assert result.beta.tolist() == [-0.2]


def test_cli(data1, tmp_path):
    path = tmp_path / "data1.tsv"
    data1.to_csv(path, sep="\t", index=False)
    script = Path(__file__).parents[1] / "examples" / "rfuwas_regulatory_prepare.py"
    proc = subprocess.run(
        [
            sys.executable,
            str(script),
            "--input",
            str(path),
            "--format",
            "tsv",
            "--release",
            "synthetic",
            "--outdir",
            str(tmp_path / "out"),
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )
    assert proc.returncode == 0, proc.stderr
    assert len(pd.read_csv(tmp_path / "out" / "rfu_qtl.tsv", sep="\t")) == 2
