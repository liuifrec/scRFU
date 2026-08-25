from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from examples.synthetic_scaling_benchmark import generate_receptors
from scrfu import pp, tl


@pytest.mark.parametrize("unique_fraction", [0.01, 0.5, 1.0])
def test_generated_reconstruction_and_deduplication_properties(unique_fraction: float) -> None:
    rows = generate_receptors(2_000, unique_fraction=unique_fraction, random_state=17)
    canonical = pp.canonicalize_receptor_table(rows)
    pp.validate_receptor_table(canonical)
    unique = canonical.drop_duplicates("cdr3aa").copy()
    unique["rfu_label"] = [f"RFU{index % 17}" for index in range(len(unique))]
    rebuilt = canonical.merge(
        unique[["cdr3aa", "rfu_label"]],
        on="cdr3aa",
        how="left",
        sort=False,
        validate="many_to_one",
    )
    assert len(rebuilt) == len(canonical)
    assert rebuilt["input_row_id"].tolist() == canonical["input_row_id"].tolist()
    assert rebuilt.groupby("cdr3aa", dropna=False)["rfu_label"].nunique().le(1).all()


def test_generated_downstream_mathematical_properties() -> None:
    rows = generate_receptors(5_000, unique_fraction=0.25, random_state=29)
    rows["rfu_label"] = [f"RFU{value % 31}" for value in pd.factorize(rows["cdr3aa"])[0]]
    rows["pass_thr"] = np.arange(len(rows)) % 5 != 0
    proportion = tl.rfu_pseudobulk(
        rows, sample_key="sample", normalize="proportion", weighting="cell"
    ).matrix
    per_1000 = tl.rfu_pseudobulk(
        rows, sample_key="sample", normalize="counts_per_1000", weighting="cell"
    ).matrix
    clr = tl.rfu_pseudobulk(rows, sample_key="sample", normalize="clr", weighting="cell").matrix
    assert proportion.sum(axis=1).to_numpy() == pytest.approx(np.ones(len(proportion)))
    assert per_1000.sum(axis=1).to_numpy() == pytest.approx(np.full(len(per_1000), 1000.0))
    assert clr.mean(axis=1).to_numpy() == pytest.approx(np.zeros(len(clr)), abs=1e-12)
    overlap = tl.rfu_overlap(proportion, metric="cosine", min_abundance=0).matrix
    assert overlap.to_numpy() == pytest.approx(overlap.to_numpy().T)
    assert np.diag(overlap) == pytest.approx(np.ones(len(overlap)))
    coupling = tl.rfu_phenotype_coupling(rows, phenotype_key="phenotype")
    for column in (
        "phenotype_specific_proportion",
        "normalized_phenotype_entropy",
        "phenotype_specificity",
        "dominant_phenotype_fraction",
    ):
        assert coupling[column].between(0, 1).all()


def test_generated_random_state_is_bitwise_deterministic() -> None:
    first = generate_receptors(1_000, random_state=101)
    second = generate_receptors(1_000, random_state=101)
    pd.testing.assert_frame_equal(first, second, check_exact=True)
