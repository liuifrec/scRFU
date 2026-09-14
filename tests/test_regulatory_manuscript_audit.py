import json

import pytest

from examples.regulatory_manuscript_audit import (
    complete_stage,
    evidence_tier,
    location,
    workbook_audit,
)


def test_early_counts_cannot_masquerade_as_completed_stage(tmp_path):
    (tmp_path / "completion.json").write_text('{"status": "running"}')
    (tmp_path / "evidence_counts.json").write_text("{}")
    with pytest.raises(FileNotFoundError):
        complete_stage(tmp_path, ["evidence_counts.json", "figure.pdf"])
    assert json.loads((tmp_path / "completion.json").read_text())["status"] == "running"
    (tmp_path / "figure.pdf").write_bytes(b"synthetic")
    complete_stage(tmp_path, ["evidence_counts.json", "figure.pdf"])
    assert json.loads((tmp_path / "completion.json").read_text())["status"] == "complete"


def test_conditional_layers_are_required_for_strongest_tier():
    assert evidence_tier(True, True, False, False).startswith("A_")
    assert evidence_tier(True, False, False, True).startswith("B_")
    assert evidence_tier(False, False, True, True).startswith("C_")
    # A nearby coloc or target annotation is deliberately not an input to this rule.


def test_strict_trb_and_flanking_coordinates():
    assert location("7", 142713241) == "TRB_locus"
    assert location("7", 142847296) == "TRB_flanking"
    assert location("6", 142713241) == "other"


def test_wrong_workbook_checksum_stops_inference(tmp_path):
    source = tmp_path / "sources/matos"
    source.mkdir(parents=True)
    (source / "supplementary_tables.xlsx").write_bytes(b"not the published workbook")
    with pytest.raises(ValueError, match="checksum"):
        workbook_audit(tmp_path, set())
