"""Independent synthetic plotting fixtures shared by preparation/drawing tests."""

import runpy
from pathlib import Path

import pytest


@pytest.fixture
def vdj_rows():
    example = Path(__file__).parents[1] / "examples/vdj_plotting.py"
    return runpy.run_path(str(example))["synthetic_receptors"]()
