import pandas as pd
import pytest

pytest.importorskip("pyarrow")

from examples.regulatory_locus_background import summarize


def test_background_counts_opportunities_and_preserves_conditional_membership(tmp_path):
    file = tmp_path / "synthetic.parquet"
    pd.DataFrame(
        {
            "variant_id": ["7:142400000[b38]A,G"] * 2 + ["7:142400010[b38]A,AT"],
            "phenotype_id": ["gene1", "gene2", "gene1"],
            "af": [0.8, 0.8, 0.1],
            "pval_nominal": [0.1, 1e-8, 0.3],
            "start_distance": [5, -10, 15],
        }
    ).to_parquet(file)
    result = summarize(file, pd.DataFrame({"variant_key": ["7:142400010:A:AT"]}))
    assert result.n_targets.tolist() == [2, 1]
    assert result.independent.tolist() == [False, True]
    assert result.variant_class.tolist() == ["SNV", "indel"]
    assert result.maf.tolist() == pytest.approx([0.2, 0.1])
