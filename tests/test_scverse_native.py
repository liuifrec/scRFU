from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData, read_h5ad

from scrfu.adapters import adapt_airr_dataframe
from scrfu.pp import canonicalize_receptor_table
from scrfu.rfu import RFURunResult
from scrfu.scverse import concat_scrfu, validate_scrfu_schema
from scrfu.tl import assign_rfu, call_rfu_table

ak = pytest.importorskip("awkward")


def _rfu_dir(path: Path) -> Path:
    path.mkdir()
    (path / "RFU.R").write_text("AssignRFUs <- function(ff) {}\n", encoding="utf-8")
    (path / "trimerMDSfit_small.Rdata").write_bytes(b"trimer")
    (path / "km5000noMax.Rdata").write_bytes(b"centers")
    return path


def _wrapper(path: Path) -> Path:
    path.write_text("#!/usr/bin/env Rscript\n", encoding="utf-8")
    return path


def _fake_run(self: object, features: pd.DataFrame, **kwargs: object) -> RFURunResult:
    del self, kwargs
    rows = features.copy().reset_index(drop=True)
    cdr3 = rows["cdr3aa"].astype("string")
    eligible = cdr3.str.startswith("C", na=False)
    unique = {value: index + 1 for index, value in enumerate(sorted(cdr3[eligible].unique()))}
    rows["eligibility_status"] = "ineligible_cdr3_not_starting_c"
    rows.loc[eligible, "eligibility_status"] = "eligible"
    rows["unique_sequence_id"] = pd.Series(pd.NA, index=rows.index, dtype="string")
    rows["rfu_id"] = pd.Series(pd.NA, index=rows.index, dtype="Int64")
    rows["rfu_label"] = pd.Series(pd.NA, index=rows.index, dtype="string")
    rows["rfu_score"] = np.nan
    rows["pass_thr"] = pd.Series(pd.NA, index=rows.index, dtype="boolean")
    for index in rows.index[eligible]:
        identifier = unique[str(cdr3[index])]
        rows.loc[index, "unique_sequence_id"] = f"sequence_{identifier:08d}"
        rows.loc[index, "rfu_id"] = identifier
        rows.loc[index, "rfu_label"] = f"RFU{identifier}"
        rows.loc[index, "rfu_score"] = 0.95 - identifier / 100
        rows.loc[index, "pass_thr"] = identifier % 2 == 1
    rows["rfu_status"] = rows["eligibility_status"]
    rows.loc[eligible & rows["pass_thr"].fillna(False), "rfu_status"] = "assigned_threshold_pass"
    rows.loc[eligible & ~rows["pass_thr"].fillna(False), "rfu_status"] = "assigned_below_threshold"
    return RFURunResult(
        rows,
        "",
        "",
        0,
        metadata={
            "original_row_count": len(rows),
            "eligible_row_count": int(eligible.sum()),
            "unique_query_count": len(unique),
            "reconstructed_output_row_count": len(rows),
            "rfu_threshold": 0.6,
        },
    )


def _chain(
    locus: str,
    junction_aa: str | None,
    v_call: str | None,
    *,
    productive: bool = True,
) -> dict[str, object]:
    return {
        "locus": locus,
        "junction_aa": junction_aa,
        "v_call": v_call,
        "productive": productive,
        "j_call": "TRBJ1-1" if locus == "TRB" else "TRAJ1",
    }


def _adata() -> AnnData:
    records = [
        [],
        [_chain("TRA", "CAVR", "TRAV1")],
        [_chain("TRB", "CASSA", "TRBV1")],
        [_chain("TRA", "CAVS", "TRAV2"), _chain("TRB", "CASSB", "TRBV2")],
        [_chain("TRB", "CASSA", "TRBV1"), _chain("TRB", "CASSC", "TRBV3")],
        [_chain("TRB", "CASSA", "TRBV9")],
        [_chain("TRB", "CASSD", None)],
        [_chain("TRB", "CASSN", "TRBV4", productive=False)],
        [_chain("TRB", "ASSBAD", "TRBV5")],
        [_chain("IGH", "CARIGH", "IGHV1")],
        [_chain("TRA", "CASSA", "TRAV8"), _chain("TRB", None, "TRBV6")],
        [{"junction_aa": "CUNKNOWN", "productive": True}],
    ]
    adata = AnnData(obs=pd.DataFrame(index=[f"c{i}" for i in range(len(records))]))
    adata.obsm["airr"] = ak.Array(records)
    return adata


@pytest.fixture
def configured(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> dict[str, Path]:
    monkeypatch.setattr("scrfu.backends.rfu_repo.RFURepoBackend.run", _fake_run)
    return {
        "rfu_dir": _rfu_dir(tmp_path / "rfu"),
        "wrapper_r_path": _wrapper(tmp_path / "wrapper.R"),
    }


def _assign(adata: AnnData, configured: dict[str, Path], **kwargs: object) -> None:
    assign_rfu(adata, **configured, **kwargs)


def test_native_chain_alignment_multichain_and_no_airr_mutation(
    configured: dict[str, Path],
) -> None:
    adata = _adata()
    before = ak.to_list(adata.obsm["airr"])
    _assign(adata, configured)

    assert ak.to_list(adata.obsm["airr"]) == before
    assert ak.to_list(ak.num(adata.obsm["scrfu"], axis=1)) == [len(row) for row in before]
    result = ak.to_list(adata.obsm["scrfu"])
    assert result[0] == []
    assert result[1][0]["assignment_status"] == "non_target_locus"
    assert result[2][0]["eligible"] is True
    assert result[4][0]["unique_sequence_id"] == result[5][0]["unique_sequence_id"]
    assert result[4][0]["rfu_id"] == result[5][0]["rfu_id"]
    assert result[7][0]["assignment_status"] == "nonproductive_chain"
    assert result[8][0]["eligible"] is False
    assert result[9][0]["assignment_status"] == "non_target_locus"
    assert result[10][1]["assignment_status"] == "missing_sequence"
    assert result[11][0]["assignment_status"] == "non_target_locus"
    assert validate_scrfu_schema(adata)["chain_count"] == sum(map(len, before))
    assert "rfu_dir" not in adata.uns["scrfu"]
    assert adata.uns["scrfu"]["runtime_paths_stored"] is False


def test_non_target_only_object_needs_no_external_rfu_configuration(tmp_path: Path) -> None:
    adata = AnnData(obs=pd.DataFrame(index=["empty", "tra-only"]))
    adata.obsm["airr"] = ak.Array([[], [_chain("TRA", "CAVR", "TRAV1")]])
    assign_rfu(adata)
    assert ak.to_list(adata.obsm["scrfu"]["assignment_status"]) == [
        [],
        ["non_target_locus"],
    ]
    assert adata.uns["scrfu"]["reference_identity"] == "unconfigured:no-productive-trb"
    path = tmp_path / "no-trb.h5ad"
    adata.write_h5ad(path)
    assert validate_scrfu_schema(read_h5ad(path))["chain_count"] == 1


def test_generic_airr_adapter_flattens_current_awkward_chains_without_expansion() -> None:
    adata = _adata()
    expected = sum(len(row) for row in ak.to_list(adata.obsm["airr"]))
    adapted = adapt_airr_dataframe(
        adata,
        productive_only=False,
        primary_chain=False,
    )
    assert len(adapted.receptors) == expected
    assert adapted.receptors["source_slot"].tolist() == [
        index for row in ak.to_list(adata.obsm["airr"]) for index in range(len(row))
    ]


def test_table_anndata_and_mudata_scientific_parity(configured: dict[str, Path]) -> None:
    mudata = pytest.importorskip("mudata")
    adata = _adata()
    explicit = assign_rfu(adata, inplace=False, **configured)
    assert explicit is not None
    canonical = canonicalize_receptor_table(
        explicit.table_result.per_row[
            [
                "input_row_id",
                "cell_id",
                "chain",
                "cdr3aa",
                "v_call",
                "productive",
                "source_adapter",
                "source_row_id",
            ]
        ]
    )
    table_result = call_rfu_table(canonical, **configured)
    columns = [
        "input_row_id",
        "eligibility_status",
        "unique_sequence_id",
        "rfu_id",
        "rfu_label",
        "rfu_score",
        "pass_thr",
    ]
    pd.testing.assert_frame_equal(
        explicit.table_result.per_row[columns].reset_index(drop=True),
        table_result.per_row[columns].reset_index(drop=True),
        check_dtype=False,
    )

    gex = AnnData(X=np.ones((1, 2)), obs=pd.DataFrame(index=["gex-only-cell"]))
    mdata = mudata.MuData({"airr": _adata(), "gex": gex})
    assign_rfu(mdata, **configured)
    assert ak.almost_equal(mdata.mod["airr"].obsm["scrfu"], explicit.chain_records)
    assert "scrfu" not in mdata.mod["gex"].obsm


def test_cell_summary_is_explicit_and_ambiguity_aware(configured: dict[str, Path]) -> None:
    adata = _adata()
    _assign(adata, configured)
    assert not any(column.startswith("scrfu_") for column in adata.obs)

    adata = _adata()
    _assign(adata, configured, cell_summary=True)
    assert adata.obs.loc["c2", "scrfu_summary_status"] == "unambiguous"
    assert adata.obs.loc["c4", "scrfu_summary_status"] == "ambiguous_multiple_rfus"
    assert pd.isna(adata.obs.loc["c4", "scrfu_rfu_id"])
    assert adata.uns["scrfu"]["summary_policy"] == "ambiguity_aware"


@pytest.mark.parametrize(
    ("policy", "status"),
    [
        ("highest_score", "selected_highest_score"),
        ("highest_threshold_score", "selected_highest_threshold_score"),
        ("first_eligible", "selected_first_eligible"),
    ],
)
def test_explicit_multi_trb_summary_policies(
    configured: dict[str, Path], policy: str, status: str
) -> None:
    adata = _adata()
    assign_rfu(adata, cell_summary=True, summary_policy=policy, **configured)
    assert adata.obs.loc["c4", "scrfu_summary_status"] == status
    assert pd.notna(adata.obs.loc["c4", "scrfu_rfu_id"])


def test_primary_summary_requires_chain_indices(configured: dict[str, Path]) -> None:
    with pytest.raises(KeyError, match="index_chains"):
        assign_rfu(
            _adata(),
            cell_summary=True,
            summary_policy="primary_vdj",
            **configured,
        )
    with pytest.raises(ValueError, match="Unknown cell-summary policy"):
        assign_rfu(_adata(), summary_policy="invented", **configured)


def test_index_chains_does_not_change_chain_assignments(configured: dict[str, Path]) -> None:
    scirpy = pytest.importorskip("scirpy")
    adata = _adata()
    before = assign_rfu(adata, inplace=False, **configured)
    scirpy.pp.index_chains(adata)
    after = assign_rfu(adata, inplace=False, **configured)
    assert ak.almost_equal(before.chain_records, after.chain_records)

    assign_rfu(
        adata,
        cell_summary=True,
        summary_policy="primary_vdj",
        **configured,
    )
    assert "scrfu_summary_status" in adata.obs


def test_anndata_h5ad_roundtrip_and_subsetting(tmp_path: Path, configured: dict[str, Path]) -> None:
    import anndata as ad

    adata = _adata()
    _assign(adata, configured, cell_summary=True)
    expected = ak.to_list(adata.obsm["scrfu"])
    path = tmp_path / "native.h5ad"
    adata.write_h5ad(path)
    loaded = ad.read_h5ad(path)
    assert ak.to_list(loaded.obsm["scrfu"]) == expected
    assert validate_scrfu_schema(loaded)["status"] == "ok"
    assert loaded.uns["scrfu"]["summary_policy"] == "ambiguity_aware"

    selected = loaded[["c2", "c4", "c8"]].copy()
    assert ak.to_list(selected.obsm["scrfu"]) == [expected[2], expected[4], expected[8]]
    assert validate_scrfu_schema(selected)["observation_count"] == 3


def test_mudata_h5mu_roundtrip_and_subsetting(tmp_path: Path, configured: dict[str, Path]) -> None:
    mudata = pytest.importorskip("mudata")
    airr = _adata()
    gex = AnnData(X=np.ones((airr.n_obs, 2)), obs=airr.obs.copy())
    mdata = mudata.MuData({"airr": airr, "gex": gex})
    assign_rfu(mdata, cell_summary=True, **configured)
    path = tmp_path / "native.h5mu"
    mdata.write_h5mu(path)
    loaded = mudata.read_h5mu(path)
    assert validate_scrfu_schema(loaded)["status"] == "ok"
    assert "scrfu" not in loaded.mod["gex"].obsm
    subset = loaded[["c2", "c4"]].copy()
    assert validate_scrfu_schema(subset)["observation_count"] == 2


def test_safe_concat_preserves_compatible_annotations_and_rejects_reference_mismatch(
    configured: dict[str, Path],
) -> None:
    first = _adata()[:5].copy()
    second = _adata()[5:].copy()
    first.obs_names = [f"a_{name}" for name in first.obs_names]
    second.obs_names = [f"b_{name}" for name in second.obs_names]
    _assign(first, configured)
    _assign(second, configured)
    combined = concat_scrfu([first, second])
    assert validate_scrfu_schema(combined)["observation_count"] == 12

    second.uns["scrfu"]["reference_identity"] = "sha256:different"
    with pytest.raises(ValueError, match="incompatible RFU references|different RFU references"):
        concat_scrfu([first, second])


def test_x_none_and_expression_matrix_remain_untouched(configured: dict[str, Path]) -> None:
    receptor_only = _adata()
    assert receptor_only.X is None
    _assign(receptor_only, configured)
    assert receptor_only.X is None

    source = _adata()
    matrix = np.arange(source.n_obs * 3).reshape(source.n_obs, 3)
    adata = AnnData(X=matrix.copy(), obs=source.obs.copy())
    adata.obsm["airr"] = source.obsm["airr"]
    matrix_identity = adata.X
    before = adata.X.copy()
    _assign(adata, configured)
    assert adata.X is matrix_identity
    np.testing.assert_array_equal(adata.X, before)


def test_current_scirpy_airr_end_to_end_roundtrip(
    tmp_path: Path, configured: dict[str, Path]
) -> None:
    scirpy = pytest.importorskip("scirpy")
    from scirpy.io import AirrCell, from_airr_cells

    cells = []
    for cell_id, cdr3 in (("cell-a", "CASSA"), ("cell-b", "CASSB")):
        cell = AirrCell(cell_id)
        chain = AirrCell.empty_chain_dict()
        chain.update(
            {
                "locus": "TRB",
                "junction": "TGTGCC",
                "junction_aa": cdr3,
                "v_call": "TRBV1",
                "j_call": "TRBJ1-1",
                "productive": True,
            }
        )
        cell.add_chain(chain)
        cells.append(cell)
    adata = from_airr_cells(cells)
    scirpy.pp.index_chains(adata)
    assign_rfu(adata, cell_summary=True, summary_policy="primary_vdj", **configured)
    path = tmp_path / "scirpy-flow.h5ad"
    adata.write_h5ad(path)
    loaded = scirpy.io.read_h5ad(path)
    assert validate_scrfu_schema(loaded)["chain_count"] == 2
    assert loaded.obs["scrfu_rfu_label"].notna().all()
