"""Small synthetic checks for the context pack's scientific boundaries."""

import h5py
import numpy as np
import pandas as pd
import pytest
from scipy.sparse import csr_matrix

from manuscript.scripts import singlecell_context_figures as pack


def test_logo_uses_unique_strings_and_modal_length_not_expanded_cell_counts():
    cells = pd.DataFrame(
        {
            "rfu": ["RFU1"] * 12,
            "cell_id": [f"c{i}" for i in range(12)],
            "cdr3aa": ["ACD", "AED"] + ["ACDD"] * 10,
            "donor": ["a", "b"] + ["c"] * 10,
            "v_call": ["V1"] * 12,
        }
    )
    seq, freq, support = pack.logo_tables(cells, ["RFU1"])
    assert support.iloc[0][
        ["cells", "unique_aa", "modal_length", "logo_unique_aa", "logo_cells"]
    ].tolist() == [12, 3, 3, 2, 2]
    assert not seq.set_index("cdr3aa").loc["ACDD", "included_in_logo"]
    middle = freq[freq.position.eq(2)].set_index("amino_acid")
    assert middle.loc["C", "frequency"] == middle.loc["E", "frequency"] == 0.5
    assert freq.groupby("position").sequence_count.sum().eq(2).all()


def test_support_ranking_does_not_select_largest_single_donor_clone():
    data = pd.DataFrame(
        {"rfu": ["huge"] * 100 + ["broad"] * 20, "donor": ["a"] * 100 + ["a"] * 10 + ["b"] * 10}
    )
    assert pack.rank_rfus(data).rfu.tolist() == ["broad", "huge"]


def matching_fixture():
    rows = []
    for donor in range(6):
        for group, count in [("RFU1", 3), ("RFU2", 6)]:
            for i in range(count):
                rows.append(
                    {
                        "cell_id": f"{donor}-{group}-{i}",
                        "rfu": group,
                        "donor": f"W{donor}",
                        "donor_id": str(donor),
                        "tissue": "blood",
                        "cell_type": "source_type",
                        "library_id": str(donor),
                        "v_call": "V1",
                    }
                )
    return pd.DataFrame(rows)


def test_matching_is_balanced_deterministic_and_never_reuses_cells():
    source = matching_fixture()
    config = {
        "case_rfu": "RFU1",
        "minimum_cells_per_group_per_donor": 3,
        "match_columns": ["donor_id", "tissue", "cell_type", "library_id", "v_call"],
    }
    matched, strata = pack.match_cells(source, config, "seed")
    shuffled, _ = pack.match_cells(source.sample(frac=1, random_state=42), config, "seed")
    assert matched.cell_id.tolist() == shuffled.cell_id.tolist()
    assert len(matched) == 36 and matched.cell_id.is_unique and matched.included_in_de.all()
    assert strata.matched_per_group.eq(3).all()
    # A different V gene cannot be silently treated as a matched background.
    source.loc[source.rfu.eq("RFU2"), "v_call"] = "V2"
    with pytest.raises(ValueError, match="Insufficient donor replication"):
        pack.match_cells(source, config, "seed")


def test_raw_reader_preserves_selected_row_and_gene_order(tmp_path):
    path = tmp_path / "source.h5ad"
    source = np.array([[1, 0, 2], [100, 100, 100], [3, 4, 0], [0, 1, 5]], dtype=float)
    matrix = csr_matrix(source)
    with h5py.File(path, "w") as f:
        x = f.create_group("raw/X")
        x.attrs["encoding-type"] = "csr_matrix"
        for key in ["data", "indices", "indptr"]:
            x.create_dataset(key, data=getattr(matrix, key))
        for location in ["raw/var", "var"]:
            v = f.create_group(location)
            v.attrs["_index"] = "_index"
            v.create_dataset("_index", data=np.array(["g1", "g2", "g3"], dtype="S"))
            v.create_dataset("feature_name", data=np.array(["A", "B", "C"], dtype="S"))
        f["var"].create_dataset("feature_is_filtered", data=[False, False, False])
    cells = pd.DataFrame(
        {
            "source_row": [3, 0, 2],
            "cell_id": ["c3", "c0", "c2"],
            "donor": ["W1", "W1", "W1"],
            "group": ["case", "case", "background"],
        }
    )
    counts, meta, genes, _ = pack.aggregate_selected_counts(path, cells)
    assert genes.gene_id.tolist() == ["g1", "g2", "g3"]
    assert counts[meta.group.eq("case")].tolist() == [[1, 1, 7]]
    assert counts[meta.group.eq("background")].tolist() == [[3, 4, 0]]
    with h5py.File(path, "a") as f:
        f["raw/X/data"][0] = 1.25
    with pytest.raises(ValueError, match="not integer counts"):
        pack.aggregate_selected_counts(path, cells)


def test_de_uses_donor_pairs_excludes_tcr_and_has_correct_direction():
    rng = np.random.default_rng(12)
    counts = rng.poisson(100, (16, 6))
    counts[::2, 0] += 150
    samples = pd.DataFrame(
        {
            "donor": np.repeat([f"W{i}" for i in range(8)], 2),
            "group": ["RFU_positive", "background"] * 8,
            "cells": [1000] * 16,
        }
    )
    genes = pd.DataFrame(
        {
            "gene_id": [f"g{i}" for i in range(6)],
            "gene": ["GENE1", "GENE2", "GENE3", "GENE4", "TRBV1", "FILTERED"],
            "source_filtered": [False] * 5 + [True],
        }
    )
    forward, _ = pack.paired_de(counts, samples, genes)
    swapped = samples.assign(
        group=samples.group.map({"RFU_positive": "background", "background": "RFU_positive"})
    )
    reverse, _ = pack.paired_de(counts, swapped, genes)
    assert forward.donors.eq(8).all()
    assert forward.loc[0, "mean_paired_log2cpm_difference"] > 0
    np.testing.assert_allclose(
        forward.mean_paired_log2cpm_difference, -reverse.mean_paired_log2cpm_difference
    )
    np.testing.assert_allclose(forward.q_value, reverse.q_value, equal_nan=True)
    assert forward.test_status.tolist()[-2:] == ["TCR_gene_excluded", "source_filtered"]
    assert forward.iloc[-2:].p_value.isna().all()
    assert forward.loc[forward.test_status.eq("tested"), "q_value"].between(0, 1).all()
