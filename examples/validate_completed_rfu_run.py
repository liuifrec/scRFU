"""Validate completed RFU outputs and manifests without invoking the RFU backend."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from scrfu.completed_run import validate_completed_rfu_run


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--evidence-manifest", type=Path)
    parser.add_argument("--verify-evidence-inputs", action="store_true")
    parser.add_argument("--expected-threshold", type=float)
    parser.add_argument("--expected-run-id")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    expected = {}
    if args.expected_threshold is not None:
        expected["rfu_threshold"] = args.expected_threshold
    if args.expected_run_id is not None:
        expected["run_id"] = args.expected_run_id
    report = validate_completed_rfu_run(
        args.run_dir,
        evidence_manifest=args.evidence_manifest,
        verify_evidence_inputs=args.verify_evidence_inputs,
        expected_provenance=expected,
    )
    text = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
