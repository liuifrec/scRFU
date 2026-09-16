import numpy as np
import pandas as pd
import pytest

from manuscript.scripts.rp1_14_repair_qc import conflict_count, stratified_sample, vector_audit


def test_vector_audit_preserves_order_values_and_missing_assignments():
    old = pd.DataFrame({"rfu": [1, 2, np.nan], "max_cor": [0.59, 0.61, np.nan]})
    fixed = old.rename(columns={"rfu": "rfu_numeric"})
    assert vector_audit(old, fixed)["rows"] == 3
    with pytest.raises(ValueError, match="altered"):
        vector_audit(old, fixed.iloc[::-1])
    with pytest.raises(ValueError, match="altered"):
        vector_audit(old, fixed.assign(max_cor=[0.60, 0.61, np.nan]))


def test_stratified_parity_is_unique_deterministic_and_result_blind():
    rows = []
    for donor in ("A", "B"):
        for visit in (1, 2, 3):
            for compartment in ("CD4", "CD8"):
                for status in ("AA_label_repaired", "AA_label_unchanged"):
                    for score in (0.4, 0.58, 0.595, 0.605, 0.63, 0.7):
                        for i in range(2):
                            rows.append(
                                dict(
                                    donor=donor,
                                    visit=visit,
                                    compartment=compartment,
                                    repair_status=status,
                                    max_cor=score,
                                    rfu=f"RFU{i + 1}",
                                    aminoAcid=f"synthetic_{len(rows)}",
                                    source_row=i,
                                    productive_canonical=True,
                                )
                            )
    frame = pd.DataFrame(rows)
    a, strata, _ = stratified_sample(frame, target=200)
    b, _, _ = stratified_sample(frame.sample(frac=1, random_state=6), target=200)
    assert len(a) == a.aminoAcid.nunique() == 200
    assert set(a.aminoAcid) == set(b.aminoAcid)
    for key in ("donor", "visit", "compartment", "repair_status", "score_bin"):
        assert a[key].nunique() == strata[key].nunique()
    assert strata.selected_unique_aa.sum() == 200


def test_repeated_receptor_conflicts_are_groups_not_observations():
    frame = pd.DataFrame({"aa": ["x", "x", "x", "y"], "rfu": [1, 2, 2, 1]})
    assert conflict_count(frame, "aa", "rfu") == 1
    assert conflict_count(frame.assign(rfu=1), "aa", "rfu") == 0
