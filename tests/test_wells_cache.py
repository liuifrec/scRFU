from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import pytest

from examples.wells_atlas_workflow import extract_wells_features, main
from scrfu.rfu import RFURunResult
from scrfu.wells import (
    UnsupportedWellsH5ADLayout,
    load_wells_receptor_cache,
    prepare_wells_receptor_cache,
    read_wells_receptors_h5ad,
)


def _source(path: Path) -> Path:
    cells = pd.Index(["c1", "c2", "c3", "c4"])
    obs = pd.DataFrame(
        {
            "donor_id": ["d1", "d1", "d2", "d2"],
            "sample_id": ["s1", "s1", "s2", "s3"],
            "unused": [1, 2, 3, 4],
        },
        index=cells,
    )
    adata = ad.AnnData(X=np.arange(20, dtype=np.float32).reshape(4, 5), obs=obs)
    adata.raw = adata.copy()
    adata.layers["counts"] = np.arange(20, dtype=np.float32).reshape(4, 5)
    adata.uns["unrelated_large_field"] = np.arange(1000)
    adata.uns["TCR_IR"] = pd.DataFrame(
        {
            "TCR-IR_VDJ_1_locus": ["TRB", "TRB", "TRA", "TRB"],
            "TCR-IR_VDJ_1_junction_aa": ["CASSA", "CASSB", "CAV", "CASSA"],
            "TCR-IR_VDJ_1_v_call": ["TRBV1", "TRBV2", "TRAV1", "TRBV1"],
            "TCR-IR_VDJ_1_productive": [True, True, True, False],
        },
        index=cells,
    )
    adata.write_h5ad(path)
    return path


def _fake_rfu_dir(path: Path) -> Path:
    path.mkdir()
    (path / "RFU.R").write_text("AssignRFUs <- function(ff) {}\n", encoding="utf-8")
    (path / "trimerMDSfit_small.Rdata").write_bytes(b"trimer")
    (path / "km5000noMax.Rdata").write_bytes(b"centers")
    return path


def _fake_run(self, features: pd.DataFrame, **kwargs: object) -> RFURunResult:
    del self, kwargs
    rows = features.copy().reset_index(drop=True)
    rows["input_row_id"] = range(len(rows))
    rows["eligibility_status"] = "eligible"
    codes, _ = pd.factorize(rows["cdr3aa"], sort=False)
    rows["unique_sequence_id"] = [f"sequence_{code:08d}" for code in codes]
    rows["rfu_id"] = codes + 1
    rows["rfu_label"] = [f"RFU{code + 1}" for code in codes]
    rows["rfu_score"] = 0.9
    rows["pass_thr"] = True
    rows["rfu_status"] = "assigned_threshold_pass"
    return RFURunResult(
        rows,
        "",
        "",
        0,
        metadata={
            "run_id": "synthetic-cache-parity",
            "chunking_enabled": True,
            "chunk_size": 2,
            "chunk_count": 1,
            "completed_chunk_count": 1,
            "reused_chunk_count": 0,
            "recomputed_chunk_count": 1,
            "invalidated_chunk_count": 0,
            "failed_chunk_count": 0,
            "total_elapsed_seconds": 0.01,
            "unique_query_count": int(rows["cdr3aa"].nunique()),
            "original_row_count": len(rows),
            "reconstructed_output_row_count": len(rows),
            "run_manifest_path": None,
        },
    )


def test_targeted_extraction_does_not_read_x_raw_layers_or_unrelated_uns(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    source = _source(tmp_path / "wells.h5ad")
    original = h5py.Dataset.__getitem__

    def guarded_getitem(dataset: h5py.Dataset, key: object) -> object:
        forbidden = ("/X", "/raw", "/layers", "/uns/unrelated_large_field")
        if dataset.name == "/X" or dataset.name.startswith(forbidden[1:]):
            raise AssertionError(f"targeted reader accessed forbidden payload {dataset.name}")
        return original(dataset, key)

    monkeypatch.setattr(h5py.Dataset, "__getitem__", guarded_getitem)
    data = read_wells_receptors_h5ad(source, obs_columns=["donor_id", "sample_id"], max_cells=3)

    assert data.atlas_shape == (4, 5)
    assert data.obs.index.tolist() == ["c1", "c2", "c3"]
    assert data.obs.columns.tolist() == ["donor_id", "sample_id"]
    assert data.tcr_ir.index.equals(data.obs.index)
    assert data.tcr_ir["TCR-IR_VDJ_1_junction_aa"].tolist() == ["CASSA", "CASSB", "CAV"]


def test_selected_obs_names_preserve_requested_order(tmp_path: Path) -> None:
    source = _source(tmp_path / "wells.h5ad")

    data = read_wells_receptors_h5ad(
        source, obs_columns=["donor_id"], selected_obs_names=["c4", "c1"]
    )

    assert data.obs.index.tolist() == ["c4", "c1"]
    assert data.tcr_ir.index.tolist() == ["c4", "c1"]


def test_cache_creation_reload_and_source_fingerprint(tmp_path: Path) -> None:
    source = _source(tmp_path / "wells.h5ad")
    cache = tmp_path / "cache"

    manifest = prepare_wells_receptor_cache(
        source, cache, obs_columns=["donor_id", "sample_id"], max_cells=3
    )
    loaded = load_wells_receptor_cache(cache, source_h5ad=source)

    assert manifest["cache_schema_version"] == 1
    assert manifest["source_atlas_dimensions"] == [4, 5]
    assert manifest["selected_metadata_columns"] == ["donor_id", "sample_id"]
    assert manifest["tcr_row_count"] == 3
    assert manifest["source_h5ad"]["fingerprint_algorithm"].startswith("sha256-")
    assert (cache / "tcr_ir.tsv.gz").is_file()
    assert (cache / "obs_metadata.tsv.gz").is_file()
    assert loaded.obs.index.tolist() == ["c1", "c2", "c3"]
    assert loaded.tcr_ir.index.equals(loaded.obs.index)


def test_changed_source_invalidates_cache(tmp_path: Path) -> None:
    source = _source(tmp_path / "wells.h5ad")
    cache = tmp_path / "cache"
    prepare_wells_receptor_cache(source, cache)
    with source.open("ab") as stream:
        stream.write(b"changed")

    with pytest.raises(ValueError, match="cache is stale"):
        load_wells_receptor_cache(cache, source_h5ad=source)


def test_unsupported_tcr_ir_hdf5_layout_fails_clearly(tmp_path: Path) -> None:
    source = _source(tmp_path / "wells.h5ad")
    with h5py.File(source, "r+") as handle:
        handle["uns/TCR_IR"].attrs["encoding-type"] = "dict"

    with pytest.raises(UnsupportedWellsH5ADLayout, match="unsupported encoding-type 'dict'"):
        read_wells_receptors_h5ad(source)


def test_h5ad_and_cache_workflow_rfu_result_parity(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    source = _source(tmp_path / "wells.h5ad")
    cache = tmp_path / "cache"
    prepare_wells_receptor_cache(source, cache, obs_columns=["donor_id"])
    rfu_dir = _fake_rfu_dir(tmp_path / "RFU")
    monkeypatch.setattr("scrfu.backends.rfu_repo.RFURepoBackend.run", _fake_run)

    common = ["--rfu-dir", str(rfu_dir), "--chunk-size", "2", "--primary-chain"]
    h5ad_out = tmp_path / "h5ad-out"
    cache_out = tmp_path / "cache-out"
    assert main(["--input", str(source), "--outdir", str(h5ad_out), *common]) == 0
    assert main(["--input", str(cache), "--outdir", str(cache_out), *common]) == 0

    for name in (
        "extracted_trb.tsv.gz",
        "unique_sequence_map.tsv.gz",
        "rfu_results_per_sequence.tsv.gz",
        "rfu_results_per_cell.tsv.gz",
    ):
        left = pd.read_csv(h5ad_out / name, sep="\t")
        right = pd.read_csv(cache_out / name, sep="\t")
        pd.testing.assert_frame_equal(left, right)
    h5ad_manifest = json.loads((h5ad_out / "run_manifest.json").read_text())
    cache_manifest = json.loads((cache_out / "run_manifest.json").read_text())
    assert h5ad_manifest["input_kind"] == "h5ad_targeted_reader"
    assert cache_manifest["input_kind"] == "prepared_wells_cache"


def test_cache_and_h5ad_adapter_feature_parity(tmp_path: Path) -> None:
    source = _source(tmp_path / "wells.h5ad")
    cache = tmp_path / "cache"
    prepare_wells_receptor_cache(source, cache)

    direct = read_wells_receptors_h5ad(source)
    prepared = load_wells_receptor_cache(cache)
    direct_features, _ = extract_wells_features(direct, primary_chain=True)
    cache_features, _ = extract_wells_features(prepared, primary_chain=True)

    pd.testing.assert_frame_equal(direct_features, cache_features)
