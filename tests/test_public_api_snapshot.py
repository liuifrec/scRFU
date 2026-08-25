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
