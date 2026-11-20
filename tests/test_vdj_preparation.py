from __future__ import annotations

import copy
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

import scrfu


@pytest.mark.parametrize(
    "name",
    [
        "prepare_receptor_plot_data",
        "prepare_rfu_embedding",
        "prepare_rfu_vj_usage",
        "prepare_rfu_cdr3_spectratype",
        "prepare_rfu_membership",
        "prepare_rfu_longitudinal_presence",
    ],
)
def test_preparation_public_imports(name):
    assert callable(getattr(scrfu.tl, name))
    assert name in scrfu.tl.__all__


def test_literal_vj_and_explicit_denominators(vdj_rows):
    result = scrfu.tl.prepare_rfu_vj_usage(vdj_rows, groupby="rfu_label")
    assert result.denominators.to_dict() == {"mock_rfu_A": 7, "mock_rfu_B": 5}
    assert result.v_counts.loc["mock_rfu_A", "TRBV1"] == 5
    assert result.j_counts.loc["mock_rfu_B", "TRBJ2-3"] == 2
    assert result.literal_pairs.value.sum() == 12
    assert result.pairs.value.sum() == 12
    assert result.v_counts.sum(axis=1).equals(result.denominators)


def test_other_categories_conserve_complete_mass_and_break_ties(vdj_rows):
    a = scrfu.tl.prepare_rfu_vj_usage(vdj_rows, top_n=1)
    b = scrfu.tl.prepare_rfu_vj_usage(vdj_rows.sample(frac=1, random_state=8), top_n=1)
    assert_frame_equal(a.v_counts, b.v_counts)
    assert a.v_counts.iloc[0].to_dict() == {"TRBV1": 5, "Other V": 7}
    assert a.j_counts.iloc[0].to_dict() == {"TRBJ1-1": 7, "Other J": 5}
    tie = vdj_rows.iloc[:2].copy()
    tie["v_call"] = ["TRBV10", "TRBV2"]
    assert list(scrfu.tl.prepare_rfu_vj_usage(tie, top_n=1).v_counts) == ["TRBV2", "Other V"]


def test_cell_receptor_read_units_differ_with_known_weights(vdj_rows):
    cell = scrfu.tl.prepare_rfu_vj_usage(vdj_rows)
    unique = scrfu.tl.prepare_rfu_vj_usage(
        vdj_rows, unit="receptor", identity_cols=["cdr3aa", "v_call"]
    )
    reads = scrfu.tl.prepare_receptor_plot_data(vdj_rows, unit="read", read_count_col="read_count")
    read = scrfu.tl.prepare_rfu_vj_usage(reads)
    assert cell.denominators.iloc[0] == 12
    assert unique.denominators.iloc[0] == 5
    assert read.denominators.iloc[0] == 34
    assert read.v_counts.iloc[0]["TRBV1"] == 19


def test_missing_calls_are_not_hidden_in_other(vdj_rows):
    vdj_rows.loc[0, "v_call"] = None
    vdj_rows.loc[1, "j_call"] = ""
    result = scrfu.tl.prepare_rfu_vj_usage(vdj_rows, top_n=1)
    assert result.v_counts.iloc[0]["Missing V"] == 1
    assert result.j_counts.iloc[0]["Missing J"] == 1
    assert result.denominators.iloc[0] == 12
    result = scrfu.tl.prepare_rfu_vj_usage(vdj_rows.drop(columns="j_call"))
    assert result.j_counts.iloc[0].to_dict() == {"Missing J": 12}
    assert result.audit["missing_columns"] == ["j_call"]


@pytest.mark.parametrize("value", [None, "", pd.NA])
def test_missing_cell_ids_rejected(vdj_rows, value):
    vdj_rows.loc[0, "cell_id"] = value
    with pytest.raises(ValueError, match="identifiers"):
        scrfu.tl.prepare_receptor_plot_data(vdj_rows)


def test_multiple_chains_require_policy_and_saved_primary(vdj_rows):
    rows = vdj_rows.iloc[:2].copy()
    rows["cell_id"] = "same_cell"
    rows["primary"] = [True, False]
    rows.loc[1, "cdr3aa"] = "CQQF"
    with pytest.raises(ValueError, match="unique"):
        scrfu.tl.prepare_receptor_plot_data(rows)
    with pytest.raises(ValueError, match="cannot count all chains"):
        scrfu.tl.prepare_receptor_plot_data(rows, multi_chain="all")
    primary = scrfu.tl.prepare_receptor_plot_data(
        rows, multi_chain="primary", primary_col="primary"
    )
    assert len(primary.rows) == 1 and primary.rows.cdr3aa.iloc[0] == "CAAF"
    all_chains = scrfu.tl.prepare_receptor_plot_data(
        rows, multi_chain="all", unit="receptor", identity_cols=["cdr3aa", "v_call"]
    )
    assert scrfu.tl.prepare_rfu_vj_usage(all_chains).denominators.iloc[0] == 2
    rows["primary"] = True
    with pytest.raises(ValueError, match="Exactly one"):
        scrfu.tl.prepare_receptor_plot_data(rows, multi_chain="primary", primary_col="primary")


def test_duplicate_read_records_are_not_summed_silently(vdj_rows):
    rows = pd.concat([vdj_rows, vdj_rows.iloc[:1]], ignore_index=True)
    with pytest.raises(ValueError, match="record IDs"):
        scrfu.tl.prepare_receptor_plot_data(
            rows, unit="read", multi_chain="all", read_count_col="read_count"
        )


@pytest.mark.parametrize("value", [-1, 0.5, np.inf, np.nan])
def test_invalid_read_counts_rejected(vdj_rows, value):
    vdj_rows["read_count"] = vdj_rows.read_count.astype(float)
    vdj_rows.loc[0, "read_count"] = value
    with pytest.raises(ValueError, match="nonnegative integers"):
        scrfu.tl.prepare_receptor_plot_data(vdj_rows, unit="read", read_count_col="read_count")


def test_bulk_reads_do_not_require_cells_but_do_require_record_ids(vdj_rows):
    rows = vdj_rows.drop(columns="cell_id")
    result = scrfu.tl.prepare_receptor_plot_data(rows, unit="read", read_count_col="read_count")
    assert scrfu.tl.prepare_rfu_vj_usage(result).denominators.iloc[0] == 34
    with pytest.raises(ValueError, match="input_row_id"):
        scrfu.tl.prepare_receptor_plot_data(
            rows.drop(columns="input_row_id"), unit="read", read_count_col="read_count"
        )


def test_locus_productivity_and_qualification_are_explicit(vdj_rows):
    vdj_rows.loc[0, "chain"] = "TRA"
    vdj_rows.loc[1, "productive"] = False
    vdj_rows.loc[2, "pass_thr"] = False
    result = scrfu.tl.prepare_receptor_plot_data(vdj_rows)
    assert len(result.rows) == 11
    result = scrfu.tl.prepare_receptor_plot_data(
        vdj_rows, productive_only=True, qualified_only=True
    )
    assert len(result.rows) == 9
    assert result.audit["input_rows"] == 12
    all_loci = scrfu.tl.prepare_receptor_plot_data(vdj_rows, chain=None)
    assert len(all_loci.rows) == 12


def test_receptor_weighting_rejects_ambiguous_gene_mapping(vdj_rows):
    vdj_rows.loc[1, "j_call"] = "TRBJ2-7"
    with pytest.raises(ValueError, match="ambiguous plotted values"):
        scrfu.tl.prepare_rfu_vj_usage(vdj_rows, unit="receptor", identity_cols=["cdr3aa", "v_call"])
    result = scrfu.tl.prepare_rfu_vj_usage(
        vdj_rows, unit="receptor", identity_cols=["junction", "v_call", "j_call"]
    )
    assert result.denominators.iloc[0] == 6


def test_identity_definitions_and_tuple_collisions(vdj_rows):
    aa = scrfu.tl.prepare_rfu_membership(vdj_rows)
    nt = scrfu.tl.prepare_rfu_membership(vdj_rows, identity_cols=["v_call", "j_call", "junction"])
    assert len(aa.receptors) == 5 and len(nt.receptors) == 6
    assert aa.receptors.abundance.sum() == nt.receptors.abundance.sum() == 12
    rows = vdj_rows.iloc[:2].copy()
    rows["x"] = ["a|b", "a"]
    rows["y"] = ["c", "b|c"]
    assert len(scrfu.tl.prepare_rfu_membership(rows, identity_cols=["x", "y"]).receptors) == 2


def test_membership_conflicts_and_mixed_references(vdj_rows):
    vdj_rows.loc[1, "rfu_label"] = "different_mock_group"
    with pytest.raises(ValueError, match="multiple RFU"):
        scrfu.tl.prepare_rfu_membership(vdj_rows)
    result = scrfu.tl.prepare_rfu_membership(vdj_rows, allow_multiple_rfus=True)
    assert len(result.edges) == 6
    vdj_rows["reference_identity"] = ["reference_a"] * 6 + ["reference_b"] * 6
    with pytest.raises(ValueError, match="incompatible reference"):
        scrfu.tl.prepare_rfu_membership(vdj_rows, allow_multiple_rfus=True)


def test_confirmed_edges_and_donor_counts(vdj_rows):
    result = scrfu.tl.prepare_rfu_membership(vdj_rows, donor_col="donor")
    lookup = result.receptors.set_index("node_id").identity
    observed = {(lookup[e.receptor_node], e.rfu) for e in result.edges.itertuples()}
    expected = set(
        zip(zip(vdj_rows.cdr3aa, vdj_rows.v_call, strict=True), vdj_rows.rfu_label, strict=True)
    )
    assert observed == expected
    assert result.rfus.abundance.sum() == 12
    assert sorted(result.receptors.donors) == [1, 2, 2, 2, 2]
    unique = scrfu.tl.prepare_rfu_membership(vdj_rows, unit="receptor")
    assert unique.receptors.abundance.sum() == 5


def test_spectratype_known_bins_and_unique_weights(vdj_rows):
    result = scrfu.tl.prepare_rfu_cdr3_spectratype(vdj_rows)
    assert result.counts.iloc[0].tolist() == [5, 4, 3]
    assert result.bin_edges.tolist() == [3.5, 4.5, 5.5, 6.5]
    unique = scrfu.tl.prepare_rfu_cdr3_spectratype(
        vdj_rows, unit="receptor", identity_cols=["cdr3aa", "v_call"]
    )
    assert unique.counts.iloc[0].tolist() == [1, 2, 2]
    binned = scrfu.tl.prepare_rfu_cdr3_spectratype(vdj_rows, bins=[3, 5, 6])
    assert binned.counts.iloc[0].tolist() == [5, 7]


@pytest.mark.parametrize("sequence", ["cassf", "CA*F", "CAXF", "CA F", 123])
def test_spectratype_invalid_sequences_are_explicit(vdj_rows, sequence):
    vdj_rows["cdr3aa"] = vdj_rows.cdr3aa.astype(object)
    vdj_rows.loc[0, "cdr3aa"] = sequence
    with pytest.raises(ValueError, match="Invalid amino-acid"):
        scrfu.tl.prepare_rfu_cdr3_spectratype(vdj_rows)
    result = scrfu.tl.prepare_rfu_cdr3_spectratype(vdj_rows, invalid="exclude")
    assert result.denominators.iloc[0] == 11
    assert result.input_denominators.iloc[0] == 12
    assert result.excluded.iloc[0]["invalid"] == 1


def test_spectratype_missing_and_truncated_bins(vdj_rows):
    vdj_rows.loc[0, "cdr3aa"] = None
    result = scrfu.tl.prepare_rfu_cdr3_spectratype(vdj_rows)
    assert result.excluded.iloc[0]["missing"] == 1
    with pytest.raises(ValueError, match="cover every"):
        scrfu.tl.prepare_rfu_cdr3_spectratype(vdj_rows, bins=[3, 5])


def test_existing_embedding_is_aligned_by_id_and_unmodified(vdj_rows):
    meta = vdj_rows.set_index("cell_id")
    coords = pd.DataFrame(np.arange(24).reshape(12, 2), index=meta.index)
    original = coords.copy()
    result = scrfu.tl.prepare_rfu_embedding(coordinates=coords, metadata=meta.iloc[::-1])
    assert result.metadata.index.tolist() == coords.index.tolist()
    assert_frame_equal(result.metadata, meta)
    assert np.array_equal(result.coordinates, coords)
    assert_frame_equal(coords, original)
    missing = meta.iloc[:-1]
    with pytest.raises(ValueError, match="missing from metadata"):
        scrfu.tl.prepare_rfu_embedding(coordinates=coords, metadata=missing)
    allowed = scrfu.tl.prepare_rfu_embedding(
        coordinates=coords, metadata=missing, allow_missing_metadata=True
    )
    assert allowed.audit["missing_metadata_cells"] == 1
    assert allowed.metadata.iloc[-1].isna().all()


@pytest.mark.parametrize(
    "kind",
    [
        "duplicate_coordinate",
        "duplicate_metadata",
        "missing_id",
        "extra_metadata",
        "nan_coordinate",
        "wrong_order",
        "no_ids",
    ],
)
def test_embedding_join_failures(vdj_rows, kind):
    meta = vdj_rows.set_index("cell_id")
    coords = pd.DataFrame(np.zeros((12, 2)), index=meta.index)
    kw = {"coordinates": coords, "metadata": meta}
    if kind == "duplicate_coordinate":
        coords.index = ["x"] * 12
    elif kind == "duplicate_metadata":
        meta.index = ["x"] * 12
    elif kind == "missing_id":
        coords.index = [None, *coords.index[1:]]
    elif kind == "extra_metadata":
        meta.index = ["unknown", *meta.index[1:]]
    elif kind == "nan_coordinate":
        coords.iloc[0, 0] = np.nan
    elif kind == "wrong_order":
        kw["cell_ids"] = coords.index[::-1]
    elif kind == "no_ids":
        kw["coordinates"] = coords.to_numpy()
    with pytest.raises(ValueError):
        scrfu.tl.prepare_rfu_embedding(**kw)


def test_anndata_embedding_never_reads_expression(vdj_rows):
    class NoExpression:
        obs = vdj_rows.set_index("cell_id")
        obsm = {"X_umap": np.arange(24).reshape(12, 2)}

        @property
        def X(self):
            raise AssertionError("Expression must not be accessed")

    result = scrfu.tl.prepare_rfu_embedding(NoExpression())
    assert result.audit["embedding_computed"] is False
    with pytest.raises(KeyError, match="no embedding is computed"):
        scrfu.tl.prepare_rfu_embedding(NoExpression(), basis="missing")


def test_longitudinal_orders_counts_and_union(vdj_rows):
    result = scrfu.tl.prepare_rfu_longitudinal_presence(
        vdj_rows,
        sample_col="visit",
        sample_order=["after", "before"],
        donor_col="donor",
        identity_cols=["cdr3aa", "v_call"],
    )
    assert result.sample_order == ("after", "before")
    assert result.denominators.to_numpy().tolist() == [[3, 3], [3, 3]]
    flow = result.presence_flows.query("donor == 'invented_donor_a'").set_index("status")
    assert flow.identities.to_dict() == {"First only": 2, "Both": 1, "Last only": 1}
    assert flow.denominator.eq(4).all()
    assert list(result.matrices["invented_donor_a"].columns) == ["after", "before"]


def test_missing_visits_differ_from_observed_zero(vdj_rows):
    rows = vdj_rows.loc[~(vdj_rows.donor.eq("invented_donor_b") & vdj_rows.visit.eq("after"))]
    kw = dict(sample_col="visit", sample_order=["before", "after"], donor_col="donor")
    missing = scrfu.tl.prepare_rfu_longitudinal_presence(rows, **kw)
    assert missing.matrices["invented_donor_b"]["after"].isna().all()
    assert not missing.observed.loc["invented_donor_b", "after"]
    design = vdj_rows[["donor", "visit"]].drop_duplicates()
    known = scrfu.tl.prepare_rfu_longitudinal_presence(rows, sample_metadata=design, **kw)
    assert known.matrices["invented_donor_b"]["after"].eq(0).all()
    assert known.denominators.loc["invented_donor_b", "after"] == 0
    assert known.observed.loc["invented_donor_b", "after"]


@pytest.mark.parametrize(
    "order", [["before"], ["before", "before"], ["before", "unknown"], ["before", None]]
)
def test_invalid_sample_orders(vdj_rows, order):
    with pytest.raises(ValueError):
        scrfu.tl.prepare_rfu_longitudinal_presence(vdj_rows, sample_col="visit", sample_order=order)


def test_longitudinal_unassigned_mass_and_read_denominators(vdj_rows):
    vdj_rows.loc[0, "rfu_label"] = None
    data = scrfu.tl.prepare_receptor_plot_data(vdj_rows, unit="read", read_count_col="read_count")
    result = scrfu.tl.prepare_rfu_longitudinal_presence(
        data, sample_col="visit", sample_order=["before", "after"], donor_col="donor"
    )
    assert result.denominators.loc["invented_donor_a"].tolist() == [7, 6]
    assert result.audit["excluded_identity_weight"][0]["weight"] == 4
    with pytest.raises(ValueError, match="missing values"):
        scrfu.tl.prepare_rfu_longitudinal_presence(
            data, sample_col="visit", sample_order=["before", "after"], missing_identity="error"
        )


def test_preparation_is_readonly_and_metadata_joins_are_exact(vdj_rows):
    original = vdj_rows.copy(deep=True)
    meta = pd.DataFrame({"category": range(12)}, index=vdj_rows.cell_id).iloc[::-1]
    data = scrfu.tl.prepare_receptor_plot_data(vdj_rows, metadata=meta)
    assert data.rows.category.tolist() == list(range(12))
    assert_frame_equal(original, vdj_rows)
    with pytest.raises(ValueError, match="absent from metadata"):
        scrfu.tl.prepare_receptor_plot_data(vdj_rows, metadata=meta.iloc[:-1])
    obs = SimpleNamespace(obs=vdj_rows.drop(columns="cell_id").set_axis(vdj_rows.cell_id), obsm={})
    assert len(scrfu.tl.prepare_receptor_plot_data(obs, airr_key=None).rows) == 12


@pytest.mark.parametrize(
    "function",
    [
        "prepare_receptor_plot_data",
        "prepare_rfu_vj_usage",
        "prepare_rfu_cdr3_spectratype",
        "prepare_rfu_membership",
    ],
)
def test_empty_inputs_raise(function):
    with pytest.raises(ValueError, match="empty"):
        getattr(scrfu.tl, function)(pd.DataFrame())


@pytest.mark.parametrize(
    "kwargs",
    [
        {"unit": "umi"},
        {"unit": "receptor"},
        {"unit": "read"},
        {"multi_chain": "guess"},
        {"multi_chain": "primary"},
    ],
)
def test_invalid_preparation_contract(vdj_rows, kwargs):
    with pytest.raises(ValueError):
        scrfu.tl.prepare_receptor_plot_data(vdj_rows, **kwargs)


def test_native_airr_and_mudata_keep_chain_assignment_alignment():
    ak = pytest.importorskip("awkward")
    from anndata import AnnData

    records = [
        [
            {
                "locus": "TRA",
                "junction_aa": "CAAQF",
                "v_call": "TRAV1",
                "j_call": "TRAJ1",
                "productive": True,
                "saved_primary": False,
            },
            {
                "locus": "TRB",
                "junction_aa": "CAAF",
                "v_call": "TRBV1",
                "j_call": "TRBJ1-1",
                "productive": True,
                "saved_primary": True,
            },
        ],
        [
            {
                "locus": "TRB",
                "junction_aa": "CGGF",
                "v_call": "TRBV2",
                "j_call": "TRBJ2-1",
                "productive": True,
                "saved_primary": True,
            },
            {
                "locus": "TRB",
                "junction_aa": "CTTF",
                "v_call": "TRBV3",
                "j_call": "TRBJ2-3",
                "productive": True,
                "saved_primary": False,
            },
        ],
    ]
    results = []
    for chains in records:
        results.append(
            [
                {
                    "eligible": r["locus"] == "TRB",
                    "rfu_id": 1,
                    "rfu_label": "mock_A" if r["locus"] == "TRB" else None,
                    "rfu_score": 0.8,
                    "threshold_pass": r["locus"] == "TRB",
                    "assignment_status": "synthetic",
                    "unique_sequence_id": str(i),
                    "locus": r["locus"],
                    "chain_index": i,
                    "reference_identity": "synthetic_reference",
                }
                for i, r in enumerate(chains)
            ]
        )
    adata = AnnData(obs=pd.DataFrame({"group": ["first", "second"]}, index=["cell_a", "cell_b"]))
    adata.obsm["airr"] = ak.Array(records)
    adata.obsm["scrfu"] = ak.Array(results)
    adata.uns["scrfu"] = {"schema_version": "1.0", "reference_identity": "synthetic_reference"}
    original = ak.to_list(adata.obsm["airr"])
    with pytest.raises(ValueError, match="unique"):
        scrfu.tl.prepare_receptor_plot_data(adata)
    result = scrfu.tl.prepare_receptor_plot_data(
        adata, multi_chain="primary", primary_col="saved_primary"
    )
    assert result.rows.cdr3aa.tolist() == ["CAAF", "CGGF"]
    assert result.rows.group.tolist() == ["first", "second"]
    assert result.rows.rfu_label.tolist() == ["mock_A", "mock_A"]
    assert ak.to_list(adata.obsm["airr"]) == original
    mudata = pytest.importorskip("mudata")
    multi = mudata.MuData({"airr": adata})
    aligned = scrfu.tl.prepare_receptor_plot_data(
        multi, multi_chain="all", unit="receptor", identity_cols=["cdr3aa", "v_call"]
    )
    assert len(aligned.rows) == 3
    damaged = copy.deepcopy(results)
    damaged[1][1]["chain_index"] = 0
    adata.obsm["scrfu"] = ak.Array(damaged)
    with pytest.raises(ValueError, match="indices/loci"):
        scrfu.tl.prepare_receptor_plot_data(
            adata, multi_chain="all", unit="receptor", identity_cols=["cdr3aa"]
        )


def test_no_matplotlib_import_for_preparation():
    # The separate core-only distribution gate also verifies a truly absent optional stack.
    import subprocess
    import sys

    run = subprocess.run(
        [
            sys.executable,
            "-c",
            "import scrfu, sys; assert 'matplotlib' not in sys.modules; assert callable(scrfu.tl.prepare_rfu_membership)",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    assert run.returncode == 0


def test_unique_weight_unassigned_audit_counts_identities(vdj_rows):
    vdj_rows.loc[vdj_rows.cdr3aa.eq("CAAF"), "rfu_label"] = None
    result = scrfu.tl.prepare_rfu_membership(vdj_rows, unit="receptor")
    assert result.audit["excluded_unassigned_weight"] == 1
    assert len(result.receptors) == 4


def test_prepared_unit_and_identity_mismatches_rejected(vdj_rows):
    prepared = scrfu.tl.prepare_receptor_plot_data(vdj_rows, unit="cell")
    with pytest.raises(ValueError, match="unit disagrees"):
        scrfu.tl.prepare_rfu_vj_usage(prepared, unit="receptor")
    with pytest.raises(ValueError, match="identity_cols disagree"):
        scrfu.tl.prepare_rfu_vj_usage(prepared, identity_cols=["cdr3aa"])


def test_literal_gene_resolution_is_never_changed(vdj_rows):
    vdj_rows.loc[0, "v_call"] = "TRBV1*01,TRBV2*02"
    vdj_rows.loc[1, "v_call"] = " TRBV1*01 "
    result = scrfu.tl.prepare_rfu_vj_usage(vdj_rows)
    assert result.v_counts.iloc[0]["TRBV1*01,TRBV2*02"] == 1
    assert result.v_counts.iloc[0][" TRBV1*01 "] == 1


def test_declared_zero_donor_and_low_detection(vdj_rows):
    design = vdj_rows[["donor", "visit"]].drop_duplicates()
    design = pd.concat(
        [design, pd.DataFrame({"donor": ["invented_empty"] * 2, "visit": ["before", "after"]})]
    )
    result = scrfu.tl.prepare_rfu_longitudinal_presence(
        vdj_rows,
        sample_col="visit",
        sample_order=["before", "after"],
        donor_col="donor",
        sample_metadata=design,
        min_count=4,
    )
    assert result.denominators.loc["invented_empty"].eq(0).all()
    assert result.observed.loc["invented_empty"].all()
    assert result.presence_flows.identities.eq(0).all()
