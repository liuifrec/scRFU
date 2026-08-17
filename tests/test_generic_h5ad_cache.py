from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import pytest

from scrfu.cli import main as cli_main
from scrfu.io import (
    UnsupportedH5ADLayout,
    migrate_wells_receptor_cache,
    read_h5ad_dataframe,
    read_h5ad_obs,
    read_receptor_cache,
    validate_receptor_cache,
    write_receptor_cache,
)
from scrfu.pp import canonicalize_receptor_table
from scrfu.wells import prepare_wells_receptor_cache


def _h5ad(path: Path) -> Path:
    cells = pd.Index(["c1", "c2", "c3"])
    obs = pd.DataFrame(
        {
            "donor": pd.Categorical(["d1", "d1", "d2"]),
            "condition": ["a", "b", "a"],
        },
        index=cells,
    )
    data = ad.AnnData(X=np.arange(12).reshape(3, 4), obs=obs)
    data.raw = data.copy()
    data.layers["counts"] = data.X.copy()
    data.uns["table"] = pd.DataFrame({"value": [1, 2, 3]}, index=cells)
    data.obsm["airr"] = pd.DataFrame(
        {
            "cell_id": cells,
            "chain": ["TRB", "TRB", "TRA"],
            "cdr3aa": ["CASSA", "CASST", "CAVR"],
            "v_call": ["TRBV1", "TRBV2", "TRAV1"],
        },
        index=cells,
    )
    data.write_h5ad(path)
    return path


def test_generic_selective_h5ad_readers_and_categorical_obs(tmp_path: Path) -> None:
    path = _h5ad(tmp_path / "input.h5ad")

    uns, uns_info = read_h5ad_dataframe(path, location="uns", key="table", return_info=True)
    obsm = read_h5ad_dataframe(path, location="obsm", key="airr")
    obs, obs_info = read_h5ad_obs(path, columns=["donor"], return_info=True)

    assert uns["value"].tolist() == [1, 2, 3]
    assert obsm.index.tolist() == ["c1", "c2", "c3"]
    assert obs["donor"].tolist() == ["d1", "d1", "d2"]
    assert isinstance(obs["donor"].dtype, pd.CategoricalDtype)
    assert uns_info.row_count == 3
    assert obs_info.index_name == "_index"


def test_selective_reader_preserves_requested_order_and_does_not_read_x(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = _h5ad(tmp_path / "input.h5ad")
    original = h5py.Dataset.__getitem__

    def guarded(dataset: h5py.Dataset, key: object) -> object:
        if dataset.name == "/X" or dataset.name.startswith(("/raw", "/layers")):
            raise AssertionError(f"read forbidden payload {dataset.name}")
        return original(dataset, key)

    monkeypatch.setattr(h5py.Dataset, "__getitem__", guarded)
    obs = read_h5ad_obs(path, columns=["condition"], selected_names=["c3", "c1"])
    assert obs.index.tolist() == ["c3", "c1"]


def test_selective_reader_missing_and_unsupported_dataframe(tmp_path: Path) -> None:
    path = _h5ad(tmp_path / "input.h5ad")
    with pytest.raises(KeyError, match="missing uns"):
        read_h5ad_dataframe(path, location="uns", key="absent")
    with h5py.File(path, "r+") as handle:
        handle["uns/table"].attrs["encoding-type"] = "dict"
    with pytest.raises(UnsupportedH5ADLayout, match="unsupported encoding-type 'dict'"):
        read_h5ad_dataframe(path, location="uns", key="table")


def test_obsm_index_mismatch_fails_explicitly(tmp_path: Path) -> None:
    path = _h5ad(tmp_path / "input.h5ad")
    with h5py.File(path, "r+") as handle:
        handle["obsm/airr/_index"][0] = "different"
    with pytest.raises(UnsupportedH5ADLayout, match="does not exactly match"):
        read_h5ad_dataframe(path, location="obsm", key="airr")


def test_prepare_receptors_cli_from_airr_h5ad_and_validate_only(tmp_path: Path) -> None:
    path = _h5ad(tmp_path / "input.h5ad")
    cache = tmp_path / "cache"
    cli_main(
        [
            "prepare-receptors",
            "--input",
            str(path),
            "--adapter",
            "scirpy_airr",
            "--outdir",
            str(cache),
            "--metadata-columns",
            "donor",
        ]
    )
    loaded = read_receptor_cache(cache)
    assert loaded.receptors["chain"].tolist() == ["TRB", "TRB"]
    assert loaded.cell_metadata["donor"].tolist() == ["d1", "d1", "d2"]

    validation_target = tmp_path / "must-not-exist"
    cli_main(
        [
            "prepare-receptors",
            "--input",
            str(path),
            "--adapter",
            "scirpy_airr",
            "--outdir",
            str(validation_target),
            "--validate-only",
        ]
    )
    assert not validation_target.exists()


def test_prepare_receptors_cli_wells_adapter_writes_canonical_cache(tmp_path: Path) -> None:
    cells = pd.Index(["c1", "c2"])
    data = ad.AnnData(obs=pd.DataFrame({"donor": ["d1", "d2"]}, index=cells))
    data.uns["TCR_IR"] = pd.DataFrame(
        {
            "TCR-IR_VDJ_1_locus": ["TRB", "TRB"],
            "TCR-IR_VDJ_1_junction_aa": ["CASSA", "CASST"],
            "TCR-IR_VDJ_1_v_call": ["TRBV1", "TRBV2"],
            "TCR-IR_VDJ_1_productive": [True, True],
        },
        index=cells,
    )
    source = tmp_path / "wells.h5ad"
    data.write_h5ad(source)
    cache = tmp_path / "wells-cache"
    cli_main(
        [
            "prepare-receptors",
            "--input",
            str(source),
            "--adapter",
            "wells_tcr_ir",
            "--outdir",
            str(cache),
            "--metadata-columns",
            "donor",
        ]
    )
    loaded = read_receptor_cache(cache)
    assert loaded.manifest["source_adapter"] == "wells_tcr_ir"
    assert loaded.receptors["source_slot"].tolist() == ["VDJ_1", "VDJ_1"]


def _receptors() -> pd.DataFrame:
    return canonicalize_receptor_table(
        pd.DataFrame(
            {
                "cell_id": ["c1", "c2", "c2"],
                "chain": ["TRB", "TRB", "TRA"],
                "cdr3aa": ["CASSA", "CASST", "CAVR"],
                "v_call": ["TRBV1", "TRBV2", "TRAV1"],
                "productive": [True, True, True],
                "source_adapter": "synthetic",
                "source_row_id": ["a", "b", "c"],
            }
        )
    )


def test_generic_cache_roundtrip_checksums_and_portable_source_state(tmp_path: Path) -> None:
    source = tmp_path / "source.txt"
    source.write_text("source", encoding="utf-8")
    metadata = pd.DataFrame({"donor": ["d1", "d2"]}, index=["c1", "c2"])
    cache = tmp_path / "cache"
    manifest = write_receptor_cache(
        cache,
        _receptors(),
        metadata,
        source_adapter="synthetic",
        source_adapter_version="1",
        source_format="table",
        source_path=source,
        selected_metadata_columns=["donor"],
    )

    loaded = read_receptor_cache(cache)
    checked = validate_receptor_cache(cache, source_path=source)
    assert manifest["cache_schema_version"] == 2
    assert loaded.receptors["input_row_id"].tolist() == [
        "row_00000000",
        "row_00000001",
        "row_00000002",
    ]
    assert loaded.cell_metadata["donor"].tolist() == ["d1", "d2"]
    assert loaded.validation["source_status"] == "unavailable"
    assert checked["source_status"] == "unchanged"


def test_generic_cache_checksum_schema_alignment_and_overwrite_protection(tmp_path: Path) -> None:
    cache = tmp_path / "cache"
    metadata = pd.DataFrame(index=["c1", "c2"])
    write_receptor_cache(
        cache,
        _receptors(),
        metadata,
        source_adapter="synthetic",
        source_adapter_version="1",
        source_format="table",
    )
    with pytest.raises(FileExistsError, match="force=True"):
        write_receptor_cache(
            cache,
            _receptors(),
            metadata,
            source_adapter="synthetic",
            source_adapter_version="1",
            source_format="table",
        )
    with (cache / "receptors.tsv.gz").open("ab") as stream:
        stream.write(b"broken")
    assert validate_receptor_cache(cache)["checksum_valid"] is False
    with pytest.raises(ValueError, match="checksum mismatch"):
        read_receptor_cache(cache)


def test_legacy_cache_detection_is_explicit(tmp_path: Path) -> None:
    cache = tmp_path / "legacy"
    cache.mkdir()
    (cache / "tcr_ir.tsv.gz").write_bytes(b"legacy")
    (cache / "preparation_manifest.json").write_text(
        json.dumps({"cache_schema_version": 1}), encoding="utf-8"
    )
    report = validate_receptor_cache(cache)
    assert report["status"] == "legacy_wells_cache"
    with pytest.raises(ValueError, match="Legacy Wells cache"):
        read_receptor_cache(cache)


def test_legacy_wells_cache_migration_is_explicit_and_non_overwriting(tmp_path: Path) -> None:
    cells = pd.Index(["c1", "c2"])
    adata = ad.AnnData(obs=pd.DataFrame(index=cells))
    adata.uns["TCR_IR"] = pd.DataFrame(
        {
            "TCR-IR_VDJ_1_locus": ["TRB", "TRB"],
            "TCR-IR_VDJ_1_junction_aa": ["CASSA", "CASST"],
            "TCR-IR_VDJ_1_v_call": ["TRBV1", "TRBV2"],
            "TCR-IR_VDJ_1_productive": [True, True],
        },
        index=cells,
    )
    source = tmp_path / "wells.h5ad"
    adata.write_h5ad(source)
    legacy = tmp_path / "legacy"
    migrated = tmp_path / "migrated"
    prepare_wells_receptor_cache(source, legacy)

    migrate_wells_receptor_cache(legacy, migrated)
    loaded = read_receptor_cache(migrated)
    assert loaded.manifest["source_format"] == "legacy_wells_cache"
    assert loaded.receptors["cdr3aa"].tolist() == ["CASSA", "CASST"]
    assert (legacy / "tcr_ir.tsv.gz").is_file()
    with pytest.raises(FileExistsError):
        migrate_wells_receptor_cache(legacy, migrated)
