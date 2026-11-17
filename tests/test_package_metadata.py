from __future__ import annotations

import re
from importlib.metadata import version
from pathlib import Path

import scrfu

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python 3.10 compatibility
    import tomli as tomllib

ROOT = Path(__file__).parents[1]


def test_package_metadata_has_no_placeholders_and_correct_urls() -> None:
    metadata = (ROOT / "pyproject.toml").read_text()
    assert "REPLACE_ME" not in metadata
    assert "example.com" not in metadata
    assert 'Homepage = "https://github.com/liuifrec/scRFU"' in metadata
    assert 'Repository = "https://github.com/liuifrec/scRFU"' in metadata
    assert 'Issues = "https://github.com/liuifrec/scRFU/issues"' in metadata


def test_package_metadata_policy_and_citation_are_consistent() -> None:
    metadata = tomllib.loads((ROOT / "pyproject.toml").read_text())
    project = metadata["project"]
    assert project["requires-python"] == ">=3.10"
    assert project["authors"] == [{"name": "Yu-Chen Liu", "email": "liu_y@rerf.or.jp"}]
    classifiers = set(project["classifiers"])
    for minor in (10, 11, 12):
        assert f"Programming Language :: Python :: 3.{minor}" in classifiers
    citation = (ROOT / "CITATION.cff").read_text()
    assert f"version: {scrfu.__version__}" in citation
    assert "orcid" not in citation.lower()
    assert not re.search(r"^doi:", citation, re.MULTILINE)
    assert 'doi: "10.1038/s41467-024-52522-z"' in citation


def test_runtime_version_matches_declared_version_source() -> None:
    version_source = (ROOT / "src" / "scrfu" / "_version.py").read_text()
    match = re.search(r'^__version__\s*=\s*["\']([^"\']+)', version_source, re.MULTILINE)
    assert match is not None
    assert match.group(1) == "0.4.1rc1"
    assert scrfu.__version__ == match.group(1)
    assert version("scrfu") == scrfu.__version__
    assert f"## {scrfu.__version__}" in (ROOT / "CHANGELOG.md").read_text()


def test_sdist_scope_excludes_scientific_and_local_evidence() -> None:
    metadata = tomllib.loads((ROOT / "pyproject.toml").read_text())
    included = metadata["tool"]["hatch"]["build"]["targets"]["sdist"]["include"]
    for prefix in ("/manuscript", "/notebooks", "/.scrfu", "/results"):
        assert not any(path.startswith(prefix) for path in included)
    for required in (
        "/src",
        "/r",
        "/tests",
        "/examples",
        "/tools",
        "/docs",
        "/LICENSE",
        "/CITATION.cff",
    ):
        assert required in included
