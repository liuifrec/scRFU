from __future__ import annotations

import importlib
import json
from pathlib import Path


def test_stable_public_api_snapshot_has_no_removals() -> None:
    snapshot_path = Path(__file__).parents[1] / "docs" / "public_api_0.4.json"
    snapshot = json.loads(snapshot_path.read_text())
    assert snapshot["schema_version"] == 1
    assert snapshot["release_line"] == "0.4"
    stable = [entry for entry in snapshot["entries"] if entry["stability"] == "stable"]
    assert stable
    for entry in stable:
        assert entry["has_docstring"], (
            f"Stable public API lacks a docstring: {entry['namespace']}.{entry['name']}"
        )
        module = importlib.import_module(entry["namespace"])
        assert hasattr(module, entry["name"]), (
            f"Stable public API removed: {entry['namespace']}.{entry['name']}"
        )


def test_public_api_snapshot_fields_and_classifications() -> None:
    snapshot_path = Path(__file__).parents[1] / "docs" / "public_api_0.4.json"
    snapshot = json.loads(snapshot_path.read_text())
    required = {
        "namespace",
        "name",
        "signature",
        "stability",
        "result_type",
        "module",
        "has_docstring",
    }
    keys: set[tuple[str, str]] = set()
    for entry in snapshot["entries"]:
        assert required == set(entry)
        assert entry["stability"] in {"stable", "experimental", "compatibility"}
        key = (entry["namespace"], entry["name"])
        assert key not in keys
        keys.add(key)
        if entry["namespace"] == "scrfu.bcr":
            assert entry["stability"] == "experimental"


def test_scverse_0_5_snapshot_extends_0_4_without_stable_removals() -> None:
    docs = Path(__file__).parents[1] / "docs"
    old = json.loads((docs / "public_api_0.4.json").read_text())
    new = json.loads((docs / "public_api_0.5.json").read_text())
    assert new["release_line"] == "0.5"
    old_stable = {
        (entry["namespace"], entry["name"])
        for entry in old["entries"]
        if entry["stability"] == "stable"
    }
    new_entries = {(entry["namespace"], entry["name"]): entry for entry in new["entries"]}
    assert old_stable.issubset(new_entries)
    for name in ("assign_rfu", "concat_scrfu", "validate_scrfu_schema"):
        assert new_entries[("scrfu.tl", name)]["stability"] == "stable"
