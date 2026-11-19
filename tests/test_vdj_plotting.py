"""Synthetic rendering contracts; no biological/private fixtures."""

from __future__ import annotations

import runpy
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from matplotlib.axes import Axes
from matplotlib.figure import Figure

import scrfu


@pytest.fixture(autouse=True)
def close_figures():
    yield
    plt.close("all")


def _umap(rows, **kwargs):
    meta = rows.set_index("cell_id")
    xy = pd.DataFrame(np.arange(24).reshape(12, 2), index=meta.index)
    return scrfu.pl.rfu_umap(coordinates=xy, metadata=meta.iloc[::-1], **kwargs)


@pytest.mark.parametrize(
    "name",
    [
        "rfu_umap",
        "rfu_vj_usage",
        "rfu_cdr3_spectratype",
        "rfu_membership_network",
        "rfu_longitudinal_presence",
    ],
)
def test_public_import_and_axes_contract(name, vdj_rows, monkeypatch):
    assert name in scrfu.pl.__all__
    assert getattr(scrfu.pl, name).__doc__

    def forbidden(*args, **kwargs):
        pytest.fail("Plotting must not save or show automatically")

    monkeypatch.setattr(plt, "show", forbidden)
    monkeypatch.setattr(Figure, "savefig", forbidden)

    def draw(ax=None):
        if name == "rfu_umap":
            return _umap(vdj_rows, ax=ax)
        kwargs = {"ax": ax}
        if name == "rfu_longitudinal_presence":
            kwargs.update(
                sample_col="visit",
                sample_order=["before", "after"],
                donor_col="donor",
                donor="invented_donor_a",
            )
        return getattr(scrfu.pl, name)(vdj_rows, **kwargs)

    assert isinstance(draw(), Axes)
    fig, axes = plt.subplots(1, 2)
    assert draw(axes[1]) is axes[1]
    assert axes[1].figure is fig
    assert not axes[0].has_data()


def test_umap_keeps_coordinates_and_display_order(vdj_rows):
    ax = _umap(
        vdj_rows, display_cells=vdj_rows.cell_id.iloc[[8, 1, 4]].tolist(), highlight=["mock_rfu_A"]
    )
    np.testing.assert_array_equal(ax.collections[0].get_offsets(), [[2, 3], [8, 9], [16, 17]])
    np.testing.assert_array_equal(ax.collections[1].get_offsets(), [[2, 3]])
    full = scrfu.tl.prepare_rfu_membership(vdj_rows)
    assert full.rfus.abundance.sum() == 12


def test_umap_qualification_does_not_infer_threshold(vdj_rows):
    vdj_rows["eligible"] = True
    vdj_rows["pass_thr"] = vdj_rows.pass_thr.astype(object)
    vdj_rows.loc[0, "pass_thr"] = None
    vdj_rows.loc[1, "pass_thr"] = False
    ax = _umap(vdj_rows, qualification_col="pass_thr", eligible_col="eligible")
    assert {t.get_text() for t in ax.get_legend().texts} == {
        "Qualified",
        "Not qualified",
        "Missing qualification",
    }
    vdj_rows.loc[2, "eligible"] = False
    with pytest.raises(ValueError, match="ineligible"):
        _umap(vdj_rows, qualification_col="pass_thr", eligible_col="eligible")


@pytest.mark.parametrize(
    "kwargs",
    [
        {"display_cells": ["absent"]},
        {"highlight": ["absent"]},
        {"palette": {}},
        {"color": "absent"},
        {"display_cells": []},
    ],
)
def test_umap_bad_display_arguments(vdj_rows, kwargs):
    with pytest.raises(ValueError):
        _umap(vdj_rows, **kwargs)


def test_vj_values_and_ribbons_keep_other_counts(vdj_rows):
    result = scrfu.tl.prepare_rfu_vj_usage(vdj_rows, groupby="rfu_label", top_n=1)
    ax = scrfu.pl.rfu_vj_usage(result, colorbar=False)
    values = np.asarray(ax.collections[0].get_array()).reshape(2, 2)
    np.testing.assert_allclose(values.sum(axis=1), 1)
    np.testing.assert_allclose(values[:, 1], [2 / 7, 1])
    ribbon = scrfu.pl.rfu_vj_usage(result, style="ribbon", group="mock_rfu_A")
    assert ribbon.get_ylim() == (0, 1)
    assert any(t.get_text() == "Other V" for t in ribbon.texts)
    assert result.denominators.sum() == 12


@pytest.mark.parametrize(
    "kwargs",
    [
        {"style": "invalid"},
        {"segment": "d"},
        {"style": "ribbon"},
        {"style": "ribbon", "group": "absent"},
    ],
)
def test_vj_invalid_modes(vdj_rows, kwargs):
    with pytest.raises(ValueError):
        scrfu.pl.rfu_vj_usage(vdj_rows, groupby="rfu_label", **kwargs)


def test_spectratype_artists_equal_prepared_bins(vdj_rows):
    result = scrfu.tl.prepare_rfu_cdr3_spectratype(vdj_rows, bins=[3.5, 4.5, 6.5])
    ax = scrfu.pl.rfu_cdr3_spectratype(result, normalize=False)
    np.testing.assert_array_equal(ax.patches[0].get_data().values, [5, 7])
    np.testing.assert_array_equal(ax.patches[0].get_data().edges, [3.5, 4.5, 6.5])


def test_network_edge_identity_layout_and_sizes(vdj_rows):
    result = scrfu.tl.prepare_rfu_membership(vdj_rows)
    ax = scrfu.pl.rfu_membership_network(result, size_scale=7, rfu_size=80)
    assert {x.get_gid() for x in ax.lines} == {
        f"{x.receptor_node}->{x.rfu_node}" for x in result.edges.itertuples()
    }
    for artist, row in zip(ax.collections, result.receptors.itertuples(), strict=False):
        assert artist.get_gid() == row.node_id
        assert artist.get_sizes()[0] == 7 * row.abundance
    second = scrfu.pl.rfu_membership_network(
        vdj_rows.sample(frac=1, random_state=1), size_scale=7, rfu_size=80
    )
    for a, b in zip(ax.collections, second.collections, strict=True):
        np.testing.assert_array_equal(a.get_offsets(), b.get_offsets())
    selected = scrfu.pl.rfu_membership_network(result, max_per_rfu=1, layout="radial")
    assert len(selected.lines) == 2
    assert result.rfus.abundance.sum() == 12


def test_network_explicit_positions_and_multiple_memberships(vdj_rows):
    result = scrfu.tl.prepare_rfu_membership(vdj_rows)
    ids = list(result.receptors.node_id) + list(result.rfus.node_id)
    pos = {key: np.asarray([i, -i], dtype=object) for i, key in enumerate(ids)}
    ax = scrfu.pl.rfu_membership_network(result, positions=pos)
    for artist in ax.collections:
        np.testing.assert_array_equal(artist.get_offsets(), [pos[artist.get_gid()]])
    chosen = result.receptors.node_id.iloc[[3, 0]].tolist()
    subset = scrfu.pl.rfu_membership_network(result, receptor_nodes=chosen, positions=pos)
    assert [a.get_gid() for a in subset.collections[:2]] == chosen
    assert len(subset.lines) == 2
    with pytest.raises(ValueError, match="positions"):
        scrfu.pl.rfu_membership_network(result, positions={})
    vdj_rows.loc[0, "rfu_label"] = "mock_rfu_B"
    result = scrfu.tl.prepare_rfu_membership(vdj_rows, allow_multiple_rfus=True)
    assert isinstance(scrfu.pl.rfu_membership_network(result), Axes)
    with pytest.raises(ValueError, match="unique membership"):
        scrfu.pl.rfu_membership_network(result, layout="radial")


@pytest.mark.parametrize(
    "kwargs",
    [
        {"layout": "distance"},
        {"max_per_rfu": 0},
        {"receptor_nodes": ["absent"]},
        {"size_scale": -1},
        {"rfu_size": float("nan")},
        {"max_per_rfu": 1, "receptor_nodes": []},
    ],
)
def test_network_bad_display_options(vdj_rows, kwargs):
    with pytest.raises(ValueError):
        scrfu.pl.rfu_membership_network(vdj_rows, **kwargs)


def _long(rows, **kwargs):
    return scrfu.tl.prepare_rfu_longitudinal_presence(
        rows,
        sample_col="visit",
        sample_order=["after", "before"],
        donor_col="donor",
        identity_cols=["cdr3aa", "v_call"],
        **kwargs,
    )


def test_longitudinal_denominators_and_display_selection(vdj_rows):
    result = _long(vdj_rows)
    ax = scrfu.pl.rfu_longitudinal_presence(
        result,
        donor="invented_donor_a",
        measure="abundance",
        normalize=True,
        max_identities=1,
        colorbar=False,
    )
    np.testing.assert_allclose(
        np.asarray(ax.collections[0].get_array()).reshape(1, 2), [[1 / 3, 2 / 3]]
    )
    assert ax.get_xticklabels()[0].get_text() == "after (n=3)"
    alluvial = scrfu.pl.rfu_longitudinal_presence(
        result, donor="invented_donor_a", style="alluvial"
    )
    assert "n=4" in alluvial.get_ylabel()
    assert {x.get_text() for x in alluvial.get_legend().texts} == {
        "First only: 2",
        "Both: 1",
        "Last only: 1",
    }


def test_longitudinal_missing_visit_stays_missing(vdj_rows):
    rows = vdj_rows.loc[~(vdj_rows.donor.eq("invented_donor_a") & vdj_rows.visit.eq("after"))]
    result = _long(rows)
    ax = scrfu.pl.rfu_longitudinal_presence(result, donor="invented_donor_a", colorbar=False)
    assert np.ma.getmaskarray(ax.collections[0].get_array()).reshape(-1, 2)[:, 0].all()
    with pytest.raises(ValueError, match="missing visits"):
        scrfu.pl.rfu_longitudinal_presence(result, donor="invented_donor_a", style="alluvial")


@pytest.mark.parametrize(
    "kwargs",
    [
        {},
        {"donor": "unknown"},
        {"donor": "invented_donor_a", "normalize": True},
        {"donor": "invented_donor_a", "style": "alluvial", "max_identities": 1},
        {"donor": "invented_donor_a", "max_identities": 0},
    ],
)
def test_longitudinal_invalid_options(vdj_rows, kwargs):
    with pytest.raises(ValueError):
        scrfu.pl.rfu_longitudinal_presence(_long(vdj_rows), **kwargs)


def test_gallery_runs_from_scratch(tmp_path):
    module = runpy.run_path(str(Path(__file__).parents[1] / "examples" / "vdj_plotting.py"))
    module["render_gallery"](tmp_path)
    for suffix in ["svg", "png", "pdf"]:
        assert (tmp_path / f"synthetic_vdj_gallery.{suffix}").stat().st_size > 1000


@pytest.mark.parametrize(
    "preparer,plotter",
    [
        ("prepare_rfu_vj_usage", "rfu_vj_usage"),
        ("prepare_rfu_cdr3_spectratype", "rfu_cdr3_spectratype"),
        ("prepare_rfu_membership", "rfu_membership_network"),
    ],
)
def test_prepared_result_rejects_unit_relabeling(vdj_rows, preparer, plotter):
    result = getattr(scrfu.tl, preparer)(vdj_rows)
    with pytest.raises(ValueError, match="unit disagrees"):
        getattr(scrfu.pl, plotter)(result, unit="read")


def test_alluvial_zero_detection_and_unit_mismatch(vdj_rows):
    result = _long(vdj_rows, min_count=100)
    with pytest.raises(ValueError, match="No detected"):
        scrfu.pl.rfu_longitudinal_presence(result, donor="invented_donor_a", style="alluvial")
    with pytest.raises(ValueError, match="unit disagrees"):
        scrfu.pl.rfu_longitudinal_presence(result, unit="read")
