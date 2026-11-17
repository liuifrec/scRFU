"""Missing-value contract across pandas object and extension string dtypes."""

import numpy as np
import pandas as pd
import pytest

from scrfu.regulatory import normalize_regulatory_variants, regulatory_triangulation


@pytest.mark.parametrize("dtype", [object, "str", "string", "string[python]", "string[pyarrow]"])
@pytest.mark.parametrize("missing", [None, np.nan, pd.NA, "", "  "])
def test_optional_text_missing_is_none_without_stringification(dtype, missing):
    if dtype == "string[pyarrow]":
        pytest.importorskip("pyarrow")
    table = pd.DataFrame(
        {
            "variant_id": ["1:42:A:G", missing],
            "rsid": [missing, "rs7"],
            "rfu_label": ["RFU1", "RFU2"],
            "source": ["test"] * 2,
            "genome_build": ["hg38"] * 2,
        }
    )
    for field in [
        "effect_allele",
        "strand",
        "release",
        "context",
        "credible_set_id",
        "locus_id",
        "url",
        "sha256",
    ]:
        table[field] = pd.Series([missing, missing], dtype=dtype)
    # pandas 2's numpy str cast stringifies sentinels before scRFU sees them;
    # use its nullable extension dtype as the equivalent of pandas 3's str dtype.
    if dtype != "str" or int(pd.__version__.split(".")[0]) >= 3:
        table["variant_id"] = table.variant_id.astype(dtype)
        table["rsid"] = table.rsid.astype(dtype)
    else:
        for field in [
            "effect_allele",
            "strand",
            "release",
            "context",
            "credible_set_id",
            "locus_id",
            "url",
            "sha256",
        ]:
            table[field] = pd.Series([missing, missing], dtype="string")
    original = table.copy(deep=True)
    result = normalize_regulatory_variants(table, layer="rfu_qtl").set_index("rfu_label")
    pd.testing.assert_frame_equal(table, original)
    assert result.loc["RFU1", "variant_key"] == "1:42:A:G"
    assert result.loc["RFU1", "rsid"] is None
    assert result.loc["RFU2", "variant_id"] is None
    assert pd.isna(result.loc["RFU2", "variant_key"])
    assert result.loc["RFU2", "identity_status"] == "unresolved_variant_identity"
    for field in [
        "effect_allele",
        "strand",
        "release",
        "context",
        "credible_set_id",
        "locus_id",
        "url",
        "sha256",
    ]:
        assert result[field].tolist() == [None, None]


@pytest.mark.parametrize("missing", [None, np.nan, pd.NA, "", "  "])
@pytest.mark.parametrize("required", ["rfu_label", "source", "genome_build"])
def test_missing_required_field_still_rejected(missing, required):
    table = pd.DataFrame(
        {
            "variant_id": ["1:42:A:G"],
            "rfu_label": ["RFU1"],
            "source": ["test"],
            "genome_build": ["hg38"],
        }
    )
    table[required] = missing
    with pytest.raises(ValueError, match=f"nonempty {required}"):
        normalize_regulatory_variants(table, layer="rfu_qtl")


def test_blank_optional_numbers_and_missing_direction():
    common = {
        "variant_id": ["1:42:A:G"],
        "source": ["test"],
        "genome_build": ["hg38"],
        "beta": [" "],
        "pvalue": [""],
        "effect_allele": [pd.NA],
    }
    rfu = pd.DataFrame({**common, "rfu_label": ["RFU1"]})
    gene = pd.DataFrame({**common, "gene": ["GENE1"]})
    normalized = normalize_regulatory_variants(rfu, layer="rfu_qtl")
    assert normalized.beta.isna().all() and normalized.pvalue.isna().all()
    result = regulatory_triangulation(rfu, eqtl=gene)
    assert result.matched_evidence.direction_status.tolist() == ["missing_effect_allele"]
    assert result.rfu_summary.max_evidence_tier.tolist() == [2]
