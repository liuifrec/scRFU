import json

import pandas as pd
import pytest

from examples import run_regulatory_realdata_pilot as pilot
from examples.run_regulatory_realdata_pilot import evidence_sets


def test_evidence_uses_corrected_tests_or_independence_not_nominal_only():
    raw = pd.DataFrame({"variant_key": ["A", "B", "C"], "pvalue": [0.9, 0.01, 1e-12]})
    independent = pd.DataFrame({"variant_key": ["A", "D"]})
    result = evidence_sets(raw, independent, 1e-6)
    assert result["tested"] == {"A", "B", "C"}
    assert result["nominal"] == {"B", "C"}
    assert result["lookup_bonferroni"] == {"C"}
    assert result["independent"] == {"A", "D"}
    assert result["evidence"] == {"A", "C", "D"}


@pytest.mark.parametrize("partial", [False, True])
def test_pilot_actual_api_metadata_annotations_and_unknown_direction(
    tmp_path, monkeypatch, partial
):
    prepared = tmp_path / "prepared"
    results = tmp_path / "results"
    prepared.mkdir()
    results.mkdir()
    base = dict(
        variant_id="7:100:A:G",
        variant_key="7:100:A:G",
        chromosome="7",
        position=100,
        ref="A",
        alt="G",
        genome_build="GRCh38",
        source="synthetic",
        release="test",
        beta=0.2,
        pvalue=1e-10,
        source_file="synthetic.tsv",
    )
    pilot.save(
        pd.DataFrame([{**base, "rfu_label": "19"}]), prepared / "rfuwas_data1_rfuqtl_grch38.tsv"
    )
    for layer, target in (("eqtl", "gene"), ("caqtl", "peak")):
        row = {**base, target: "synthetic_target", "independent": not partial, "rank": 1, "se": 0.1}
        nominal_name = f"matos_{layer}_at_rfuqtl_variants.tsv"
        independent_name = f"matos_{layer}_independent_chr6_chr7.tsv"
        manifest_name = f"matos_{layer}_archive_manifest.json"
        independent = pd.DataFrame([row])
        if partial and layer == "caqtl":
            nominal_name = "matos_caqtl_chr7_partial_at_rfuqtl_variants.tsv"
            independent_name = "matos_caqtl_chr6_partial_independent.tsv"
            manifest_name = "matos_caqtl_partial_manifest.json"
            row["independent"] = pd.NA
            independent = independent.iloc[:0].drop(columns="source_file")
        pilot.save(pd.DataFrame([row]), prepared / nominal_name)
        pilot.save(independent, prepared / independent_name)
        (prepared / manifest_name).write_text("{}")
    (results / "rfuwas_context_provenance.json").write_text("{}")
    pilot.save(
        pd.DataFrame([dict(rfu_label="19", enrichment="CD4-enriched", Enriched_Group="TN")]),
        results / "rfuqtl_cellstate_context.tsv",
    )
    pilot.save(
        pd.DataFrame(
            [
                dict(
                    rfu_label="19",
                    description="synthetic disease",
                    phecode="001.2",
                    pvalue=1e-9,
                    effect=0.3,
                    group="autoimmune",
                )
            ]
        ),
        prepared / "rfuwas_data4_annotations.tsv",
    )
    monkeypatch.setattr(pilot, "plot_coverage", lambda *args, **kwargs: None)
    pilot.run(tmp_path, partial_caqtl=partial)
    out = results / "partial_caqtl_checkpoint" if partial else results
    hits = pilot.read(out / "regulatory_hits.tsv")
    assert set(hits.evidence_layer) == {"eqtl", "caqtl"}
    assert hits.direction_concordant.isna().all()
    ranked = pilot.read(out / "ranked_regulatory_candidates.tsv")
    assert ranked.rfu_label.tolist() == ["19"]
    assert ranked.strongest_disease_phecode.iloc[0] == "001.2"
    counts = json.loads((out / "overlap_counts.json").read_text())
    assert counts["both_evidence"] == {"variants": 1, "rfus": 1}
    if partial:
        assert counts["caqtl"]["independent"]["variants"] is None
        assert counts["caqtl"]["independent_variant_target_pairs"] is None
        assert ranked.matos_caqtl_independent.isna().all()
