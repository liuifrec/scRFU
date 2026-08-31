from __future__ import annotations

import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parents[1] / "src"))

from scrfu import __version__

project = "scRFU"
author = "Yu-Chen Liu"
copyright = "2026, Yu-Chen Liu"
version = __version__
release = __version__

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "sphinx.ext.napoleon",
]
if os.environ.get("READTHEDOCS") == "True" or os.environ.get("SCRFU_DOCS_INTERSPHINX") == "1":
    intersphinx_mapping = {
        "anndata": ("https://anndata.readthedocs.io/en/stable/", None),
        "scanpy": ("https://scanpy.readthedocs.io/en/stable/", None),
        "mudata": ("https://mudata.readthedocs.io/en/stable/", None),
        "scirpy": ("https://scirpy.scverse.org/en/stable/", None),
    }
else:
    # Offline local builds remain deterministic; hosted builds resolve the
    # authoritative inventories above.
    intersphinx_mapping = {}
autosummary_generate = True
autodoc_typehints = "description"
napoleon_numpy_docstring = True
napoleon_google_docstring = False
nitpicky = True
nitpick_ignore_regex = [
    ("py:class", r".*DataFrame.*"),
    ("py:class", r".*Series.*"),
    ("py:class", r".*Path.*"),
    ("py:class", r".*Sequence.*"),
    ("py:class", r".*Mapping.*"),
    ("py:class", r".*Any.*"),
    ("py:class", r".*Literal.*"),
    ("py:class", r".*Axes.*"),
    ("py:class", r"argparse\.ArgumentParser"),
    ("py:class", r"collections\.abc\.Callable"),
    ("py:class", r"anndata\..*\.AnnData"),
    ("py:class", r"numpy\.ndarray"),
]
suppress_warnings = [
    "toc.not_included",  # Historical planning records are intentionally outside release navigation.
    "myst.xref_missing",  # A legacy acquisition-directory label is not an API link.
]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
html_theme = "alabaster"
html_static_path: list[str] = []
