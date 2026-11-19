"""Run from an empty directory against an installed wheel or sdist, never editable."""

from __future__ import annotations

import argparse
import importlib
import importlib.util
import json
import os
import shutil
import subprocess
import sys
from importlib.metadata import version
from importlib.resources import files
from pathlib import Path
from unittest.mock import patch

import pandas as pd

import scrfu
from scrfu.doctor import doctor_report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--require-core-only", action="store_true")
    args = parser.parse_args()
    source = args.source_root.resolve()
    installed = Path(scrfu.__file__).resolve()
    assert installed.is_relative_to(Path(sys.prefix).resolve()), installed
    assert not installed.is_relative_to(source), installed
    assert not Path.cwd().resolve().is_relative_to(source)
    assert not os.environ.get("PYTHONPATH"), "Unset PYTHONPATH for this gate"
    assert scrfu.__version__ == version("scrfu")
    for module in [scrfu, scrfu.pp, scrfu.tl, scrfu.pl, scrfu.io]:
        for name in module.__all__:
            assert getattr(module, name) is not None, (module, name)
    for name in ["vdjdb", "benchmark", "comparators", "diagnostics", "longitudinal", "transfer"]:
        importlib.import_module(f"scrfu.{name}")
    if args.require_core_only:
        for name in ["matplotlib", "awkward", "scirpy", "scanpy", "mudata", "openpyxl", "pyarrow"]:
            assert importlib.util.find_spec(name) is None, f"Unexpected optional dependency: {name}"
        simple = pd.DataFrame(
            {
                "cell_id": ["invented_cell"],
                "chain": ["TRB"],
                "cdr3aa": ["CAAF"],
                "v_call": ["TRBV1"],
                "rfu_label": ["mock_A"],
            }
        )
        result = scrfu.tl.prepare_rfu_membership(simple)
        assert len(result.edges) == 1
        assert "matplotlib" not in sys.modules
        try:
            scrfu.pl.rfu_membership_network(result)
        except ImportError as exc:
            assert "plotting" in str(exc)
        else:
            raise AssertionError("Plotting must require the optional extra")
    cli = str(Path(sys.prefix) / ("Scripts/scrfu.exe" if os.name == "nt" else "bin/scrfu"))
    commands = [[cli, "--version"], [cli, "--help"]]
    commands += [
        [cli, subcommand, "--help"]
        for subcommand in [
            "doctor",
            "call-rfu",
            "prepare-wells",
            "prepare-receptors",
            "migrate-receptor-cache",
        ]
    ]
    for command in commands:
        subprocess.run(command, check=True, capture_output=True, text=True)
    for fixture in ["tutorial_receptors.tsv", "tutorial_bcr.tsv", "tutorial_vdjdb.tsv"]:
        assert files("scrfu").joinpath("data", fixture).read_bytes()
    subprocess.run(
        [sys.executable, "-m", "scrfu.tutorial", "--outdir", "synthetic"],
        check=True,
        capture_output=True,
        text=True,
    )
    manifest = json.loads(Path("synthetic/run_manifest.json").read_text())
    assert manifest["backend"] == "mock" and not manifest["rfu_provenance"]["external_rfu"]
    assert manifest["outputs"]["assigned_receptors.tsv"]["rows"] == 12
    assert manifest["outputs"]["rfu_pseudobulk_proportions.tsv"]["rows"] == 6
    assert all(item["synthetic"] for item in manifest["fixtures"].values())

    external = Path("mock-external").resolve()
    external.mkdir()
    for name in ["RFU.R", "km5000noMax.Rdata", "trimerMDSfit_small.Rdata"]:
        (external / name).write_text("AssignRFUs <- function(ff) {}" if name == "RFU.R" else "MOCK")
    with files("scrfu").joinpath("data", "tutorial_receptors.tsv").open("r") as handle:
        receptors = pd.read_csv(handle, sep="\t")
    selected = []

    def mock_run(command, **kwargs):
        wrapper = Path(command[1])
        assert wrapper.is_file()
        selected.append(wrapper)
        queries = pd.read_csv(command[command.index("--in") + 1], sep="\t")
        pd.DataFrame(
            {
                "unique_sequence_id": queries.unique_sequence_id,
                "rfu_id": 1,
                "rfu_label": "RFU1",
                "rfu_score": 0.8,
                "pass_thr": True,
                "upstream_n_miss": 0,
            }
        ).to_csv(command[command.index("--out") + 1], sep="\t", index=False)
        return subprocess.CompletedProcess(command, 0, "MOCK subprocess; not RFU validation", "")

    default = doctor_report(rfu_dir=external, verbose=True)
    assert default["wrapper"]["source"] == "package"
    override = Path("override.R").resolve()
    override.write_text("# explicitly selected mock wrapper")
    with patch("scrfu.backends.rfu_repo.subprocess.run", mock_run):
        for name, wrapper in [("default", None), ("override", override)]:
            result = scrfu.tl.call_rfu_table(
                receptors,
                rfu_dir=external,
                wrapper_r_path=wrapper,
                chunk_size=2,
                max_workers=2,
                executor="thread",
                workdir=name,
            )
            assert result.per_row.input_row_id.tolist() == receptors.input_row_id.tolist()
            expected = Path(default["wrapper"]["path"]) if wrapper is None else wrapper
            assert all(path == expected for path in selected)
            selected.clear()
    diagnosed = json.loads(
        subprocess.check_output(
            [
                cli,
                "doctor",
                "--verbose",
                "--rfu-dir",
                str(external),
                "--wrapper-r-path",
                str(override),
            ],
            text=True,
        )
    )
    assert diagnosed["wrapper"]["path"] == str(override)
    assert diagnosed["wrapper"]["source"] == "explicit"
    for missing in [None, float("nan"), pd.NA, "", " "]:
        raw = pd.DataFrame(
            {
                "variant_id": ["1:42:A:G", missing],
                "rsid": [missing, "rs7"],
                "rfu_label": ["RFU1", "RFU2"],
                "source": ["synthetic"] * 2,
                "genome_build": ["hg38"] * 2,
            }
        )
        normalized = scrfu.tl.normalize_regulatory_variants(raw, layer="rfu_qtl").set_index(
            "rfu_label"
        )
        assert normalized.loc["RFU1", "rsid"] is None
        assert normalized.loc["RFU2", "variant_id"] is None
        assert normalized.loc["RFU2", "identity_status"] == "unresolved_variant_identity"
    # Copy only the standalone CLI, not src or the repository, into the external test directory.
    shutil.copy2(source / "examples/rfuwas_regulatory_prepare.py", "regulatory_prepare.py")
    pd.DataFrame(
        {"SNP": ["7_12345_a_g"], "RFU": [1], "beta": [-0.2], "t.stat": [-2.0], "p.value": [0.03]}
    ).to_csv("data1.tsv", sep="\t", index=False)
    subprocess.run(
        [
            sys.executable,
            "regulatory_prepare.py",
            "--input",
            "data1.tsv",
            "--format",
            "tsv",
            "--release",
            "synthetic-installed-check",
            "--outdir",
            "regulatory",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    exported = pd.read_csv("regulatory/rfu_qtl.tsv", sep="\t")
    assert exported.variant_key.tolist() == ["7:12345:A:G"]
    assert exported.effect_allele.isna().all() and exported.beta.tolist() == [-0.2]
    report = {
        "version": scrfu.__version__,
        "python": sys.version,
        "pandas": pd.__version__,
        "package_file": str(installed),
        "core_only": args.require_core_only,
        "mock_wrapper_checks": "passed (not real RFU validation)",
        "status": "passed",
    }
    Path("installed-verification.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
