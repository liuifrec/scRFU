import pandas as pd
import pytest

from examples.rfuwas_regulatory_context import reference_qc


def variants():
    return pd.DataFrame(
        dict(
            variant_key=["1:1:A:G", "1:2:CG:C", "1:4:A:G", "2:1:A:G"],
            chromosome=["1", "1", "1", "2"],
            position=[1, 2, 4, 1],
            ref=["A", "CG", "A", "A"],
            genome_build=["GRCh38"] * 4,
        )
    )


def test_reference_positions_indels_mismatch_and_missing():
    result = reference_qc(
        variants(), [dict(genome="hg38", chrom="chr1", start=0, end=4, dna="acgt")]
    )
    assert result.reference_sequence.tolist() == ["A", "CG", "T", None]
    assert result.reference_checked.tolist() == [True, True, True, False]
    assert result.reference_match.iloc[:3].tolist() == [True, True, False]
    assert pd.isna(result.reference_match.iloc[3])


def test_reference_rejects_wrong_build_truncation_and_conflict():
    ref = dict(genome="hg38", chrom="chr1", start=0, end=4, dna="acgt")
    for refs in ([{**ref, "genome": "hg19"}], [{**ref, "end": 5}], [ref, {**ref, "dna": "tttt"}]):
        with pytest.raises(ValueError):
            reference_qc(variants(), refs)
    with pytest.raises(ValueError, match="GRCh38"):
        reference_qc(variants().assign(genome_build="GRCh37"), [ref])
