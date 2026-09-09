import pandas as pd
import pytest

from examples.matos_regulatory_lookup import significance_flags, source_variant_ids


def test_source_identity_keeps_alleles_and_indels():
    assert source_variant_ids({"7:123:AT:A"}) == {"7:123[b38]AT,A", "chr7:123[b38]AT,A"}


def test_tested_nominal_and_source_significance_are_distinct():
    nominal = pd.DataFrame(
        {"phenotype_id": ["A", "A", "B", "C"], "pval_nominal": [0.2, 0.001, 0.001, 0.001]}
    )
    top = pd.DataFrame(
        {
            "phenotype_id": ["A", "B", "C"],
            "qval": [0.05, 0.1, 0.01],
            "pval_nominal_threshold": [0.002, 0.002, 0.001],
            "pval_beta": [0.01] * 3,
            "pval_perm": [0.01] * 3,
        }
    )
    result = significance_flags(nominal, top)
    assert result.tested.all()
    assert result.nominal_p_lt_005.tolist() == [False, True, True, True]
    assert result.source_significant.tolist() == [False, True, False, False]
    with pytest.raises(ValueError, match="Missing phenotype threshold"):
        significance_flags(nominal, top.iloc[:2])
    with pytest.raises(ValueError, match="Ambiguous"):
        significance_flags(nominal, pd.concat([top, top]))


def test_absent_source_thresholds_are_unknown_not_negative():
    nominal = pd.DataFrame({"phenotype_id": ["A"], "pval_nominal": [1e-20]})
    result = significance_flags(nominal, None)
    assert result.tested.all()
    assert result.nominal_p_lt_005.all()
    assert result.source_significant.isna().all()
