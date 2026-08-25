from __future__ import annotations

import json
from pathlib import Path


def test_evidence_index_schema_and_portability() -> None:
    path = Path(__file__).parents[1] / "docs" / "evidence_index.json"
    index = json.loads(path.read_text())
    assert index["schema_version"] == 1
    assert isinstance(index["entries"], list)
    required = {
        "dataset_label",
        "analysis_type",
        "external_manifest_filename",
        "manifest_sha256",
        "public_source_accession",
        "result_status",
    }
    for entry in index["entries"]:
        assert required <= set(entry)
        assert not Path(entry["external_manifest_filename"]).is_absolute()
        assert len(entry["manifest_sha256"]) == 64
        int(entry["manifest_sha256"], 16)
        assert entry["result_status"] in {
            "sealed",
            "previously_verified_source_unavailable",
        }
