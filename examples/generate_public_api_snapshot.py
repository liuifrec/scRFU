"""Generate a machine-readable scRFU public API snapshot."""

from __future__ import annotations

import argparse
import importlib
import inspect
import json
from pathlib import Path

EXPERIMENTAL_TL = {
    "CohortHarmonizationResult",
    "ComparatorRepresentation",
    "FrozenRFUReference",
    "HeldOutValidationManifest",
    "LongitudinalCompartmentResult",
    "LongitudinalDesign",
    "LongitudinalDynamicsResult",
    "LongitudinalResamplingResult",
    "RFULongitudinalResult",
    "StabilityBenchmarkResult",
    "TransferCohortResult",
    "benchmark_representation_stability",
    "bootstrap_longitudinal_statistic",
    "create_heldout_validation_manifest",
    "deterministic_subsample",
    "donor_leave_one_out",
    "donor_retrieval",
    "harmonize_cohort_metadata",
    "list_comparators",
    "longitudinal_compartment_comparison",
    "longitudinal_similarity",
    "multinomial_abundance_resample",
    "permute_longitudinal_labels",
    "reference_coverage",
    "register_comparator",
    "repertoire_representation",
    "rfu_longitudinal_dynamics",
    "rfu_longitudinal_matrix",
    "shuffle_input_order",
    "summarize_longitudinal_similarity",
    "threshold_sensitivity",
    "transfer_cohort",
    "validate_frozen_reference",
    "validate_longitudinal_design",
}
COMPATIBILITY_TL = {"aggregate_rfu", "rfu_summary"}


def _result_type(obj: object) -> str | None:
    try:
        annotation = inspect.signature(obj).return_annotation
    except (TypeError, ValueError):
        return None
    if annotation is inspect.Signature.empty:
        return None
    return str(annotation)


def _entry(namespace: str, name: str, stability: str) -> dict[str, object]:
    module = importlib.import_module(namespace)
    obj = getattr(module, name)
    try:
        signature = str(inspect.signature(obj))
    except (TypeError, ValueError):
        signature = None
    return {
        "namespace": namespace,
        "name": name,
        "signature": signature,
        "stability": stability,
        "result_type": _result_type(obj),
        "module": getattr(obj, "__module__", namespace),
        "has_docstring": bool(inspect.getdoc(obj)),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("docs/public_api_0.4.json"))
    parser.add_argument("--release-line", default="0.4")
    args = parser.parse_args()
    entries: list[dict[str, object]] = []
    for namespace in ("scrfu.pp", "scrfu.adapters", "scrfu.io", "scrfu.tl", "scrfu.pl"):
        module = importlib.import_module(namespace)
        for name in sorted(module.__all__):
            obj = getattr(module, name)
            if namespace == "scrfu.adapters" and name == "ReceptorAdapter":
                continue
            if not (inspect.isfunction(obj) or inspect.isclass(obj)):
                continue
            stability = "stable"
            if namespace == "scrfu.tl" and name in EXPERIMENTAL_TL:
                stability = "experimental"
            elif namespace == "scrfu.tl" and name in COMPATIBILITY_TL:
                stability = "compatibility"
            entries.append(_entry(namespace, name, stability))
    bcr = importlib.import_module("scrfu.bcr")
    for name in sorted(bcr.__all__):
        entries.append(_entry("scrfu.bcr", name, "experimental"))
    entries.extend(
        [
            _entry("scrfu.cli", "build_parser", "stable"),
            _entry("scrfu.cli", "entrypoint", "stable"),
            _entry("scrfu.cli", "main", "stable"),
            _entry("scrfu.doctor", "doctor_report", "stable"),
        ]
    )
    payload = {
        "schema_version": 1,
        "release_line": args.release_line,
        "removal_policy": "Stable entries may not be removed without an explicit API review.",
        "entries": sorted(entries, key=lambda item: (item["namespace"], item["name"])),
    }
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
