from __future__ import annotations

from types import SimpleNamespace

import pandas as pd
import pytest

import scrfu


def _bcr_source() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "barcode": ["c1", "c1", "c1", "c2", "c2", "c3"],
            "contig_id": ["h1-low", "h1", "k1", "h2", "l2", "h3"],
            "locus": ["IGH", "heavy", "kappa", "IGH", "lambda", "IGH"],
            "junction_aa": ["CARAAA", "CARBBB", "CQQAA", "CARCCC", "CQQCC", "CARDDD"],
            "junction": ["TGT"] * 6,
            "v_gene": ["IGHV1", "IGHV1", "IGKV1", "IGHV2", "IGLV2", "IGHV3"],
            "j_gene": ["IGHJ1", "IGHJ1", "IGKJ1", "IGHJ2", "IGLJ2", "IGHJ3"],
            "constant": ["IGHM*01", "IGHG1*01", "IGKC", "IgA1", "IGLC", pd.NA],
            "is_productive": [False, True, True, True, True, pd.NA],
            "umi_count": [100, 10, 20, 4, 5, 1],
            "clone_id": ["cl1", "cl1", "cl1", "cl2", "cl2", pd.NA],
            "family_id": ["f1", "f1", "f1", "f2", "f2", pd.NA],
            "v_identity": ["99%", "0.95", "1", "90%", "bad", pd.NA],
        }
    )


def test_bcr_preparation_normalizes_and_retains_selection_provenance() -> None:
    result = scrfu.bcr.prepare_bcr_table(_bcr_source(), source_label="synthetic")
    receptors = result.receptors
    assert receptors["chain"].tolist() == ["IGH", "IGH", "IGK", "IGH", "IGL", "IGH"]
    assert receptors["isotype"].tolist()[:4] == ["IgM", "IgG1", pd.NA, "IgA1"]
    assert receptors["sequence_id"].is_unique
    assert len(receptors) == len(_bcr_source())
    selected_c1 = receptors.loc[
        receptors["cell_id"].eq("c1") & receptors["selected_for_pair"], "sequence_id"
    ]
    assert set(selected_c1) == {"h1", "k1"}
    assert not receptors.loc[receptors["sequence_id"].eq("h1-low"), "selected_for_pair"].item()
    assert result.pairs.set_index("cell_id").loc["c1", "pair_status"] == "paired"
    assert result.pairs.set_index("cell_id").loc["c3", "pair_status"] == "heavy_only"
    assert result.provenance["field_provenance"]["mutation_frequency"]["status"] == "derived"
    assert not result.provenance["tcr_rfu_assignment_permitted"]


def test_bcr_state_features_are_conservative_and_bounded() -> None:
    result = scrfu.bcr.prepare_bcr_table(_bcr_source())
    features = scrfu.bcr.bcr_state_features(result.receptors, result.pairs)
    switched = features.set_index("sequence_id")["class_switched"]
    assert switched["h1-low"] == False  # noqa: E712
    assert switched["h1"] == True  # noqa: E712
    assert pd.isna(switched["k1"])
    assert pd.isna(switched["h3"])
    assert features.loc[features["clonal_family_id"].eq("f1"), "clonal_family_size"].eq(3).all()
    diversity = features["within_family_cdr3_diversity"].dropna()
    assert diversity.between(0, 1).all()
    assert result.qc["pair_status_counts"] == {"paired": 2, "heavy_only": 1}


def test_bcr_schema_rejects_missing_or_duplicate_identifiers_and_never_assigns_tcr_rfu() -> None:
    with pytest.raises(ValueError, match="missing required"):
        scrfu.bcr.prepare_bcr_table(pd.DataFrame({"cell_id": ["c1"]}))
    duplicate = _bcr_source().assign(contig_id="same")
    with pytest.raises(ValueError, match="must be unique"):
        scrfu.bcr.prepare_bcr_table(duplicate)
    assert not hasattr(scrfu.bcr, "call_rfu")


@pytest.mark.parametrize(
    ("raw", "expected"),
    [
        ("IGHM*01", "IgM"),
        ("IgD", "IgD"),
        ("IGHG4", "IgG4"),
        ("IgA2", "IgA2"),
        ("IGHE", "IgE"),
        ("unknown", pd.NA),
    ],
)
def test_isotype_normalization(raw, expected) -> None:
    observed = scrfu.bcr.normalize_isotype(raw)
    assert pd.isna(observed) if pd.isna(expected) else observed == expected


def test_airr_anndata_and_mudata_like_bcr_routing_retains_heavy_and_light() -> None:
    airr = _bcr_source().rename(
        columns={
            "barcode": "cell_id",
            "contig_id": "sequence_id",
            "locus": "chain",
            "junction_aa": "cdr3aa",
            "v_gene": "v_call",
            "is_productive": "productive",
        }
    )
    obs = pd.DataFrame({"phenotype": ["x", "y", "z"]}, index=["c1", "c2", "c3"])
    adata = SimpleNamespace(obsm={"airr": airr}, obs=obs)
    direct = scrfu.adapters.prepare_receptors(
        adata,
        adapter="generic_airr_dataframe",
        chain=None,
        primary_chain=False,
        productive_only=False,
        metadata_columns=["phenotype"],
    )
    assert direct.receptors["chain"].tolist() == ["IGH", "IGH", "IGK", "IGH", "IGL", "IGH"]
    assert len(direct.receptors) == len(airr)
    mudata = SimpleNamespace(mod={"airr_mod": adata, "obs_mod": SimpleNamespace(obs=obs)})
    routed = scrfu.adapters.prepare_receptors(
        mudata,
        adapter="generic_airr_dataframe",
        modality="airr_mod",
        metadata_modality="obs_mod",
        chain=None,
        primary_chain=False,
        productive_only=False,
        metadata_columns=["phenotype"],
    )
    assert routed.receptors["source_row_id"].tolist() == direct.receptors["source_row_id"].tolist()
    prepared = scrfu.bcr.prepare_bcr_table(routed.receptors)
    assert prepared.pairs.set_index("cell_id").loc["c2", "light_chain"] == "IGL"
    assert routed.provenance["mudata_modality"] == "airr_mod"
