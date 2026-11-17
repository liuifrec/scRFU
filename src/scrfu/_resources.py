"""Resolve scRFU-owned resources without retaining temporary extraction paths."""

from __future__ import annotations

import hashlib
from collections.abc import Iterator
from contextlib import contextmanager
from importlib import resources
from pathlib import Path


def _wrapper_resource():
    return resources.files("scrfu").joinpath("resources", "run_rfu_repo.R")


@contextmanager
def wrapper_path(override: str | Path | None = None) -> Iterator[Path]:
    """Keep the default resource alive for the entire caller's subprocess lifetime.

    Workers resolve independently; no temporary path or open context is pickled.
    Explicit paths retain their usual filesystem semantics.
    """
    if override is not None:
        path = Path(override).expanduser().resolve()
        if not path.is_file():
            raise FileNotFoundError(f"Wrapper R script not found: {path}")
        yield path
    else:
        with resources.as_file(_wrapper_resource()) as path:
            yield path


def wrapper_details(override: str | Path | None = None) -> dict[str, str]:
    """Describe the selected wrapper; archive resources use a persistent locator."""
    resource = (
        Path(override).expanduser().resolve() if override is not None else _wrapper_resource()
    )
    with wrapper_path(override) as path:
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
    return {
        "path": str(resource)
        if isinstance(resource, Path)
        else "package:scrfu/resources/run_rfu_repo.R",
        "source": "explicit" if override is not None else "package",
        "sha256": digest,
    }
