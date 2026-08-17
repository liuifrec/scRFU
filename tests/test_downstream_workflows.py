from __future__ import annotations

import json

import pandas as pd
import pytest

from examples.benchmark_summary import main as benchmark_main
from examples.benchmark_summary import summarize_manifest
from examples.rfu_downstream_analysis import build_parser, run_analysis


def _results() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "cell_id": ["c1", "c2", "c3", "c4"],
            "cdr3aa": ["CASSA", "CASSA", "CASSB", "CASSC"],
            "v_call": ["TRBV1", "TRBV1", "TRBV2", "TRBV3"],
            "productive": [True] * 4,
            "rfu_label": ["RFU1", "RFU1", "RFU2", "RFU2"],
            "pass_thr": [True, False, True, True],
        }
    )


def test_downstream_parser_help() -> None:
    with pytest.raises(SystemExit) as error:
        build_parser().parse_args(["--help"])
    assert error.value.code == 0


def test_downstream_synthetic_end_to_end(tmp_path) -> None:
    results = tmp_path / "results.tsv.gz"
    metadata = tmp_path / "metadata.tsv.gz"
    outdir = tmp_path / "downstream"
    _results().to_csv(results, sep="\t", index=False)
    pd.DataFrame(
        {
            "cell_id": ["c1", "c2", "c3", "c4"],
            "sample": ["s1", "s1", "s2", "s2"],
            "phenotype": ["A", "B", "A", "B"],
        }
    ).to_csv(metadata, sep="\t", index=False)
    manifest = run_analysis(
        rfu_results=results,
        metadata=metadata,
        sample_key="sample",
        phenotype_key="phenotype",
        outdir=outdir,
    )
    expected = {
        "repertoire_metrics.tsv",
        "rfu_metrics.tsv",
        "rfu_pseudobulk_counts.tsv",
        "rfu_pseudobulk_proportions.tsv",
        "rfu_overlap_jaccard.tsv",
        "rfu_overlap_cosine.tsv",
        "rfu_phenotype_coupling.tsv",
        "rfu_metric_heatmap.png",
        "rfu_overlap_heatmap.png",
        "rfu_convergence.png",
        "rfu_phenotype_heatmap.png",
        "run_manifest.json",
    }
    assert expected == {path.name for path in outdir.iterdir()}
    assert manifest["join_qc"]["unmatched_result_cells"] == 0
    assert manifest["dimensions"]["pseudobulk_samples"] == 2


def test_downstream_metadata_mismatch_and_duplicate_join(tmp_path) -> None:
    results = tmp_path / "results.tsv"
    metadata = tmp_path / "metadata.tsv"
    _results().to_csv(results, sep="\t", index=False)
    pd.DataFrame({"cell_id": ["c1", "c2"], "sample": ["s1", "s1"], "phenotype": ["A", "B"]}).to_csv(
        metadata, sep="\t", index=False
    )
    with pytest.raises(ValueError, match="missing requested labels"):
        run_analysis(
            rfu_results=results,
            metadata=metadata,
            sample_key="sample",
            phenotype_key="phenotype",
            outdir=tmp_path / "out",
        )


def test_benchmark_summary_normalization_and_cli(tmp_path) -> None:
    manifest = {
        "source_atlas_dimensions": [100, 50],
        "receptor_row_count": 80,
        "eligible_row_count": 60,
        "unique_query_count": 40,
        "assigned_row_count": 55,
        "threshold_pass_row_count": 30,
        "chunk_count": 4,
        "reused_chunk_count": 2,
        "total_elapsed_seconds": 12.5,
        "max_rss_bytes": 1234,
        "source_adapter": "synthetic",
        "cache_schema_version": 2,
        "mode": "standard",
        "artifact_hashes": {"RFU.R": "abc"},
    }
    row = summarize_manifest(manifest, "dataset-a")
    assert row["input_cells"] == 100
    assert row["deduplication_ratio"] == 2
    assert row["cache_reuse_count"] == 2
    source = tmp_path / "manifest.json"
    output = tmp_path / "summary.tsv"
    source.write_text(json.dumps(manifest), encoding="utf-8")
    benchmark_main(
        [
            "--manifest",
            str(source),
            "--dataset-name",
            "dataset-a",
            "--output",
            str(output),
        ]
    )
    assert pd.read_csv(output, sep="\t").loc[0, "dataset_name"] == "dataset-a"
