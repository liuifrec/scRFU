"""Small, read-only installation diagnostics used by ``scrfu doctor``."""

from __future__ import annotations

import importlib.util
import os
import platform
import shutil
import sys
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from ._version import __version__
from .backends.rfu_repo import RFUCapabilities, RFURepoPaths
from .io import file_sha256


def _display_path(path: Path, *, verbose: bool) -> str:
    return str(path) if verbose else f"…/{path.name}"


def _writable(path: Path) -> bool:
    candidate = path.expanduser()
    while not candidate.exists() and candidate != candidate.parent:
        candidate = candidate.parent
    return candidate.is_dir() and os.access(candidate, os.W_OK)


def doctor_report(
    *,
    verbose: bool = False,
    output_dir: str | Path = ".scrfu",
    environ: Mapping[str, str] | None = None,
) -> dict[str, Any]:
    """Return redacted, read-only installation and external-backend diagnostics.

    Parameters
    ----------
    verbose
        Include full runtime paths. The default exposes only basenames.
    output_dir
        Prospective cache/output directory checked for writability. It is not
        created by this function.
    environ
        Optional environment mapping, primarily for controlled diagnostics and
        tests. Explicit arguments still take precedence in analysis APIs.
    """
    environment = os.environ if environ is None else environ
    rscript = shutil.which("Rscript")
    report: dict[str, Any] = {
        "scrfu_version": __version__,
        "python_version": platform.python_version(),
        "python_executable": _display_path(Path(sys.executable), verbose=verbose),
        "operating_system": platform.platform(),
        "rscript_available": rscript is not None,
        "rscript": _display_path(Path(rscript), verbose=verbose) if rscript else None,
        "rfu_dir_configured": bool(environment.get("RFU_DIR")),
        "rfu_dir": None,
        "rfu_required_files_present": False,
        "rfu_capability_mode": "unconfigured",
        "rfu_artifact_hashes": {},
        "plotting_available": importlib.util.find_spec("matplotlib") is not None,
        "mudata_available": importlib.util.find_spec("mudata") is not None,
        "experimental_bcr_preprocessing_available": importlib.util.find_spec("scrfu.bcr")
        is not None,
        "output_directory": _display_path(Path(output_dir), verbose=verbose),
        "output_writable": _writable(Path(output_dir)),
        "vdjdb_path_configured": bool(environment.get("VDJDB_PATH")),
        "vdjdb_release": environment.get("VDJDB_RELEASE"),
    }
    if not environment.get("RFU_DIR"):
        return report
    directory = Path(environment["RFU_DIR"]).expanduser().resolve()
    report["rfu_dir"] = _display_path(directory, verbose=verbose)
    try:
        paths = RFURepoPaths.resolve(environ=environment)
    except (FileNotFoundError, ValueError) as error:
        report["rfu_error"] = str(error) if verbose else error.args[0].splitlines()[0]
        report["rfu_capability_mode"] = "invalid"
        return report
    capabilities = RFUCapabilities.from_rfu_r(paths.rfu_r)
    mode = "standard"
    if capabilities.assign_rfus_with_map:
        mode = "standard+map_aware"
    elif not capabilities.assign_rfus:
        mode = "unsupported"
    report.update(
        {
            "rfu_required_files_present": True,
            "rfu_capability_mode": mode,
            "rfu_artifact_hashes": {
                "RFU.R": file_sha256(paths.rfu_r),
                "trimerMDSfit_small.Rdata": file_sha256(paths.trimer_rdata),
                "km5000noMax.Rdata": file_sha256(paths.km5000_rdata),
            },
        }
    )
    return report


__all__ = ["doctor_report"]
