from __future__ import annotations

import hashlib

import pandas as pd
import pytest

from scrfu.tl import (
    annotate_vdjdb,
    load_vdjdb_reference,
    normalize_vdjdb_cdr3,
    normalize_vdjdb_v_gene,
    summarize_vdjdb_evidence,
)


def _reference_frame() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "CDR3": ["CASSAAA", "CASSAAA", "CASSDDD", "CAVCCC"],
            "V": ["TRBV7-9*01", "TRBV7-9*02", "TRBV5-1*01", "TRAV1-2*01"],
            "Gene": ["TRB", "TRB", "TRB", "TRA"],
            "Epitope": ["PEPTIDE_A", "PEPTIDE_B", "PEPTIDE_A", "PEPTIDE_C"],
            "MHC A": ["HLA-A*01", "HLA-A*02", "HLA-B*07", "HLA-A*02"],
            "Score": [3, 2, 1, 3],
            "Reference": ["PMID1", "PMID2", "PMID3", "PMID4"],
        }
    )


def _queries() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "input_row_id": ["r1", "r2", "r3", "r4"],
            "unique_sequence_id": ["s1", "s1", "s2", "s3"],
            "cell_id": ["c1", "c2", "c3", "c4"],
            "cdr3aa": ["CASSAAA", "CASSAAA", "CASSDDD", "CAVCCC"],
            "v_call": ["TRBV7-9*03", "TRBV7-9*03", "TRBV5-1*01", "TRAV1-2*01"],
            "chain": ["TRB", "TRB", "TRB", "TRA"],
            "rfu_label": ["RFU1", "RFU1", "RFU1", "RFU2"],
            "pass_thr": [True, True, False, True],
        }
    )


def test_reference_minimal_alternate_columns_and_deterministic_ids() -> None:
    reference = load_vdjdb_reference(_reference_frame(), release_label="synthetic-2026-01")
    assert len(reference.table) == 4
    assert reference.table["reference_row_id"].tolist() == [
        "vdjdb_000000000",
        "vdjdb_000000001",
        "vdjdb_000000002",
        "vdjdb_000000003",
    ]
    assert reference.validation["available_chains"] == ["TRA", "TRB"]
    assert reference.provenance["release_label"] == "synthetic-2026-01"
    assert reference.provenance["unique_cdr3_count"] == 3


def test_reference_gzip_checksum_and_mismatch(tmp_path) -> None:
    path = tmp_path / "reference.tsv.gz"
    _reference_frame().to_csv(path, sep="\t", index=False, compression="gzip")
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    reference = load_vdjdb_reference(
        path,
        release_label="snapshot-1",
        expected_sha256=digest,
    )
    assert reference.provenance["sha256"] == digest
    assert reference.provenance["source_filename"] == "reference.tsv.gz"
    with pytest.raises(ValueError, match="SHA256 mismatch"):
        load_vdjdb_reference(path, release_label="snapshot-1", expected_sha256="0" * 64)


def test_reference_contract_errors_and_non_strict_validation() -> None:
    with pytest.raises(ValueError, match="release_label"):
        load_vdjdb_reference(_reference_frame(), release_label="")
    with pytest.raises(ValueError, match="required CDR3"):
        load_vdjdb_reference(pd.DataFrame({"V": ["TRBV1"]}), release_label="x")
    invalid = pd.DataFrame({"CDR3": ["CASS?"], "Epitope": ["A"]})
    non_strict = load_vdjdb_reference(invalid, release_label="x", strict=False)
    assert non_strict.validation["status"] == "error"
    valid_query = _queries().iloc[[0]]
    assert annotate_vdjdb(valid_query, non_strict, match_mode="cdr3").empty
    with pytest.raises(ValueError, match="invalid amino-acid"):
        load_vdjdb_reference(invalid, release_label="x", strict=True)


def test_reference_explicit_column_overrides_and_duplicate_rows() -> None:
    raw = pd.DataFrame(
        {"aa": ["CASSAAA", "CASSAAA"], "antigen": ["A", "A"], "locus": ["TRB", "TRB"]}
    )
    reference = load_vdjdb_reference(
        raw,
        release_label="x",
        column_mappings={"cdr3aa": "aa", "epitope": "antigen"},
    )
    assert reference.validation["duplicated_database_row_count"] == 2
    assert reference.table["reference_row_id"].is_unique


def test_normalization_is_conservative() -> None:
    assert normalize_vdjdb_cdr3(" cassaaa ") == "CASSAAA"
    assert normalize_vdjdb_v_gene("trbv7-9*01", mode="exact") == "TRBV7-9*01"
    assert normalize_vdjdb_v_gene("trbv7-9*01", mode="strip_allele") == "TRBV7-9"
    assert pd.isna(normalize_vdjdb_v_gene(pd.NA))
    with pytest.raises(ValueError, match="invalid amino-acid"):
        normalize_vdjdb_cdr3("CASS-X")


def test_exact_matching_tiers_ambiguity_and_safe_reconstruction() -> None:
    reference = load_vdjdb_reference(_reference_frame(), release_label="snapshot-1")
    evidence = annotate_vdjdb(_queries(), reference, match_mode="cdr3_v", chain="TRB")
    assert evidence["input_row_id"].tolist() == ["r1", "r1", "r2", "r2", "r3"]
    assert evidence["evidence_tier"].unique().tolist() == ["cdr3_v_exact"]
    assert evidence.loc[evidence["unique_sequence_id"].eq("s1"), "epitope"].nunique() == 2
    summaries = summarize_vdjdb_evidence(_queries(), evidence)
    assert len(summaries.row_summary) == len(_queries())
    s1 = summaries.sequence_summary.set_index("unique_sequence_id").loc["s1"]
    assert s1["antigen_ambiguity"]
    assert pd.isna(s1["dominant_epitope"])
    assert not summaries.sequence_summary.set_index("unique_sequence_id").loc[
        "s3", "has_vdjdb_evidence"
    ]


def test_exact_allele_mode_chain_mismatch_and_no_match() -> None:
    reference = load_vdjdb_reference(_reference_frame(), release_label="snapshot-1")
    exact = annotate_vdjdb(_queries(), reference, match_mode="cdr3_v", v_gene_mode="exact")
    assert set(exact["unique_sequence_id"]) == {"s2"}
    tra_query = _queries().loc[_queries()["unique_sequence_id"].eq("s3")]
    assert annotate_vdjdb(tra_query, reference, match_mode="cdr3", chain="TRB").empty
    tra = annotate_vdjdb(tra_query, reference, match_mode="cdr3", chain=None)
    assert tra["matched_chain"].tolist() == ["TRA"]
    no_match = _queries().iloc[[0]].assign(cdr3aa="CASSGGG")
    assert annotate_vdjdb(no_match, reference, match_mode="cdr3", chain="TRB").empty


def test_paired_exact_requires_real_pair_fields() -> None:
    reference = load_vdjdb_reference(_reference_frame(), release_label="snapshot-1")
    with pytest.raises(ValueError, match="genuine paired_cdr3aa"):
        annotate_vdjdb(_queries(), reference, match_mode="paired_exact")


def test_paired_exact_preserves_explicit_pair_evidence() -> None:
    reference = load_vdjdb_reference(
        pd.DataFrame(
            {
                "cdr3aa": ["CASSAAA"],
                "v_call": ["TRBV1*01"],
                "chain": ["TRB"],
                "paired_cdr3aa": ["CAVDDD"],
                "paired_v_call": ["TRAV1*01"],
                "pair_id": ["pair-reference-1"],
                "epitope": ["A"],
            }
        ),
        release_label="paired-synthetic",
    )
    query = pd.DataFrame(
        {
            "input_row_id": ["r1"],
            "unique_sequence_id": ["s1"],
            "cdr3aa": ["CASSAAA"],
            "v_call": ["TRBV1*02"],
            "chain": ["TRB"],
            "paired_cdr3aa": ["CAVDDD"],
            "paired_v_call": ["TRAV1*02"],
        }
    )
    evidence = annotate_vdjdb(query, reference, match_mode="paired_exact")
    assert evidence.loc[0, "evidence_tier"] == "paired_exact"
    assert evidence.loc[0, "query_paired_cdr3aa"] == "CAVDDD"
    assert evidence.loc[0, "matched_paired_cdr3aa"] == "CAVDDD"
    assert evidence.loc[0, "reference_paired_receptor_id"] == "pair-reference-1"
