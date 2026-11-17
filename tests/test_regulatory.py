from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

import scrfu

TARGETS = {"rfu_qtl": "rfu_label", "eqtl": "gene", "caqtl": "peak", "gwas": "trait"}


def evidence(layer="rfu_qtl", **changes):
    row = dict(
        chromosome="chr1",
        position=12345,
        ref="a",
        alt="g",
        genome_build="hg38",
        source="synthetic",
        release="v1",
        beta=0.4,
        effect_allele="g",
    )
    row[TARGETS[layer]] = {
        "rfu_qtl": "RFU1",
        "eqtl": "GENE1",
        "caqtl": "chr1:12300-12400",
        "gwas": "trait1",
    }[layer]
    row.update(changes)
    return pd.DataFrame([row])


@pytest.mark.parametrize(
    ("chrom", "expected"),
    [
        ("chr1", "1"),
        ("CHR01", "1"),
        (1.0, "1"),
        ("chrx", "X"),
        ("chrM", "MT"),
        ("MT", "MT"),
        ("GL000207.1", "GL000207.1"),
    ],
)
def test_normalization(chrom, expected):
    raw = evidence(chromosome=chrom, position="12345.0", variant_id="rs42")
    before = raw.copy(deep=True)
    result = scrfu.tl.normalize_regulatory_variants(raw, layer="rfu_qtl")
    assert result.loc[0, "variant_key"] == f"{expected}:12345:A:G"
    assert result.loc[0, "position"] == 12345
    assert result.loc[0, "rsid"] == "rs42"
    assert result.loc[0, "genome_build"] == "GRCh38"
    assert_frame_equal(raw, before)


def test_coordinate_id_and_optional_schema():
    raw = pd.DataFrame(
        [dict(variant_id="chr1:12345:a:g", rfu_label="RFU1", source="test", genome_build="GRCh38")]
    )
    result = scrfu.tl.regulatory_triangulation(raw)
    assert result.harmonized_rfu_qtl.loc[0, "variant_key"] == "1:12345:A:G"
    assert pd.isna(result.harmonized_rfu_qtl.loc[0, "beta"])
    assert result.rfu_summary.loc[0, "max_evidence_tier"] == 1
    for layer in TARGETS:
        schema = scrfu.tl.regulatory_evidence_schema(layer)
        assert TARGETS[layer] in schema.required
        assert "pip" in schema.optional


@pytest.mark.parametrize("layer", ["eqtl", "caqtl", "gwas"])
def test_single_external_layer(layer):
    result = scrfu.tl.regulatory_triangulation(evidence(), **{layer: evidence(layer)})
    assert len(result.matched_evidence) == 1
    hit = result.matched_evidence.iloc[0]
    assert hit.shared_variant_exact and hit.allele_harmonized and hit.direction_concordant
    assert result.harmonized_rfu_qtl.loc[0, f"has_{layer}"]
    assert result.rfu_summary.loc[0, "max_evidence_tier"] == (1 if layer == "gwas" else 2)


@pytest.mark.parametrize(
    ("layers", "tier", "count"),
    [
        ([], 1, 1),
        (["gwas"], 1, 2),
        (["eqtl"], 2, 2),
        (["caqtl", "gwas"], 2, 3),
        (["eqtl", "caqtl"], 3, 3),
        (["eqtl", "caqtl", "gwas"], 4, 4),
    ],
)
def test_tiers(layers, tier, count):
    result = scrfu.tl.regulatory_triangulation(
        evidence(), **{layer: evidence(layer, pip=0.96) for layer in layers}
    )
    row = result.harmonized_rfu_qtl.iloc[0]
    assert row.evidence_tier == tier and row.evidence_count == count
    assert row.eqtl_finemapped == ("eqtl" in layers)
    assert row.caqtl_finemapped == ("caqtl" in layers)
    assert not result.provenance["formal_colocalization"]
    json.dumps(result.provenance, allow_nan=False)


def test_reversal_alignment_and_disabling():
    other = evidence("eqtl", ref="G", alt="A", effect_allele="A", beta=-0.2)
    result = scrfu.tl.regulatory_triangulation(evidence(), eqtl=other)
    hit = result.matched_evidence.iloc[0]
    assert hit.alleles_reversed and not hit.shared_variant_exact
    assert hit.direction_concordant and hit.aligned_evidence_beta == 0.2
    assert hit.evidence_beta == -0.2 and hit.evidence_ref == "G"
    assert not result.harmonized_rfu_qtl.loc[0, "shared_variant_exact"]
    assert scrfu.tl.regulatory_triangulation(
        evidence(), eqtl=other, allow_allele_reversal=False
    ).matched_evidence.empty


@pytest.mark.parametrize(
    ("left", "right", "status", "harmonized"),
    [
        ({}, {"effect_allele": None}, "missing_effect_allele", False),
        ({}, {"strand": "-"}, "unsupported_negative_strand", False),
        (
            {"alt": "T", "effect_allele": "T"},
            {"alt": "T", "effect_allele": "T"},
            "palindromic_strand_unknown",
            False,
        ),
        ({}, {"beta": None}, "missing_effect_size", True),
        ({}, {"beta": 0}, "zero_effect", True),
    ],
)
def test_unknown_direction(left, right, status, harmonized):
    hit = scrfu.tl.regulatory_triangulation(
        evidence(**left), eqtl=evidence("eqtl", **right)
    ).matched_evidence.iloc[0]
    assert hit.direction_status == status
    assert hit.allele_harmonized == harmonized
    assert pd.isna(hit.direction_concordant)


def test_palindrome_explicit_forward_and_odds_ratio():
    rfu = evidence(alt="T", effect_allele="T", strand="+")
    eqtl = evidence("eqtl", ref="T", alt="A", effect_allele="A", strand="+", beta=-1)
    assert scrfu.tl.regulatory_triangulation(rfu, eqtl=eqtl).matched_evidence.loc[
        0, "direction_concordant"
    ]
    hit = scrfu.tl.regulatory_triangulation(
        evidence(), gwas=evidence("gwas", beta=None, odds_ratio=0.5)
    ).matched_evidence.iloc[0]
    assert not hit.direction_concordant
    assert hit.aligned_evidence_beta == pytest.approx(np.log(0.5))


def test_no_strand_complement_or_position_only_matching():
    external = pd.concat(
        [
            evidence("eqtl", ref="T", alt="C", effect_allele="C"),
            evidence("eqtl", alt="C", effect_allele="C"),
        ]
    )
    assert scrfu.tl.regulatory_triangulation(evidence(), eqtl=external).matched_evidence.empty


@pytest.mark.parametrize("bad", ["hg19", "GRCh37", "GRCm38"])
def test_incompatible_builds(bad):
    with pytest.raises(ValueError, match="Incompatible genome builds"):
        scrfu.tl.regulatory_triangulation(evidence(), eqtl=evidence("eqtl", genome_build=bad))
    with pytest.raises(ValueError, match="Incompatible genome builds"):
        scrfu.tl.credible_set_overlap(evidence(), evidence("eqtl", genome_build=bad))


@pytest.mark.parametrize(
    ("changes", "message"),
    [
        ({"position": 1.2}, "position"),
        ({"position": 0}, "position"),
        ({"position": True}, "position"),
        ({"position": "bad"}, "position"),
        ({"position": np.inf}, "position"),
        ({"ref": "N"}, "ref"),
        ({"alt": "A,G"}, "alt"),
        ({"alt": "A"}, "differ"),
        ({"pvalue": -0.1}, "pvalue"),
        ({"pip": 1.1}, "pip"),
        ({"se": 0}, "se"),
        ({"beta": "bad"}, "beta"),
        ({"beta": np.inf}, "beta"),
        ({"effect_allele": "C"}, "effect_allele"),
        ({"genome_build": None}, "genome_build"),
        ({"source": " "}, "source"),
        ({"rfu_label": None}, "rfu_label"),
        ({"variant_id": "1:12346:A:G"}, "conflicts"),
        ({"rsid": "bad"}, "rsid"),
        ({"strand": "forward"}, "strand"),
        ({"allele_frequency": 1.1}, "allele_frequency"),
    ],
)
def test_invalid_values(changes, message):
    with pytest.raises(ValueError, match=message):
        scrfu.tl.regulatory_triangulation(evidence(**changes))


def test_invalid_schemas_and_metadata():
    with pytest.raises(TypeError, match="DataFrame"):
        scrfu.tl.normalize_regulatory_variants([], layer="eqtl")
    with pytest.raises(ValueError, match="gene"):
        scrfu.tl.normalize_regulatory_variants(evidence(), layer="eqtl")
    with pytest.raises(ValueError, match="Variant identity"):
        scrfu.tl.normalize_regulatory_variants(pd.DataFrame([dict(gene="G")]), layer="eqtl")
    with pytest.raises(ValueError, match="Unknown evidence layer"):
        scrfu.tl.regulatory_evidence_schema("bad")
    with pytest.raises(ValueError, match="Duplicate column"):
        scrfu.tl.normalize_regulatory_variants(
            pd.DataFrame([[1, 2]], columns=["gene", "gene"]), layer="eqtl"
        )
    with pytest.raises(ValueError, match="conflicts"):
        scrfu.tl.normalize_regulatory_variants(evidence(), layer="rfu_qtl", genome_build="hg19")
    with pytest.raises(ValueError, match="input_metadata"):
        scrfu.tl.regulatory_triangulation(evidence(), input_metadata={"wrong": {}})
    for threshold in (-1, 2, np.nan):
        with pytest.raises(ValueError, match="pip_threshold"):
            scrfu.tl.regulatory_triangulation(evidence(), pip_threshold=threshold)


def test_metadata_defaults():
    raw = evidence().drop(columns=["source", "release", "genome_build"])
    meta = {
        "rfu_qtl": dict(source="user", release="v2", genome_build="hg38", url="https://example.org")
    }
    result = scrfu.tl.regulatory_triangulation(raw, input_metadata=meta)
    assert result.harmonized_rfu_qtl.loc[0, "source"] == "user"
    assert result.provenance["inputs"]["rfu_qtl"]["metadata"] == meta["rfu_qtl"]


def test_duplicates_and_determinism():
    rfu = pd.concat(
        [evidence(), evidence(), evidence(rfu_label="RFU2", position=99)], ignore_index=True
    )
    eqtl = pd.concat(
        [evidence("eqtl"), evidence("eqtl"), evidence("eqtl", gene="GENE2", beta=-0.2)],
        ignore_index=True,
    )
    first = scrfu.tl.regulatory_triangulation(rfu, eqtl=eqtl)
    second = scrfu.tl.regulatory_triangulation(
        rfu.sample(frac=1, random_state=9), eqtl=eqtl.iloc[::-1]
    )
    for field in (
        "harmonized_rfu_qtl",
        "matched_evidence",
        "unmatched_variants",
        "rfu_summary",
        "variant_summary",
        "credible_set_overlaps",
    ):
        assert_frame_equal(getattr(first, field), getattr(second, field))
    assert first.provenance == second.provenance
    assert len(first.matched_evidence) == 2
    assert first.unmatched_variants.reason.eq("duplicate_record").sum() == 2
    assert first.rfu_summary.set_index("rfu_label").loc["RFU1", "genes"] == '["GENE1", "GENE2"]'


def test_summary_does_not_construct_chain_across_variants():
    rfu = pd.concat([evidence(), evidence(position=99)])
    result = scrfu.tl.regulatory_triangulation(
        rfu,
        eqtl=evidence("eqtl"),
        caqtl=evidence("caqtl", position=99),
        gwas=evidence("gwas", position=99),
    )
    summary = result.rfu_summary.iloc[0]
    assert summary.has_eqtl and summary.has_caqtl and summary.has_gwas
    assert summary.max_evidence_tier == 2
    assert len(result.matched_evidence) == 3


def test_unmatched_and_rsid_only():
    rfu = pd.DataFrame([dict(variant_id="rs123", rfu_label="R", genome_build="hg38", source="s")])
    result = scrfu.tl.regulatory_triangulation(rfu, eqtl=evidence("eqtl", variant_id="rs123"))
    assert result.matched_evidence.empty
    assert set(result.unmatched_variants.reason) == {
        "unresolved_variant_identity",
        "no_matching_rfu_qtl",
    }
    assert result.variant_summary.empty
    assert result.rfu_summary.loc[0, "n_variants"] == 0


@pytest.mark.parametrize("external", [False, True])
def test_empty(external):
    result = scrfu.tl.regulatory_triangulation(
        pd.DataFrame(), eqtl=evidence("eqtl") if external else pd.DataFrame()
    )
    assert result.harmonized_rfu_qtl.empty and result.matched_evidence.empty
    assert result.rfu_summary.empty and result.variant_summary.empty
    assert "direction_concordant" in result.matched_evidence
    assert len(result.unmatched_variants) == int(external)
    json.dumps(result.provenance, allow_nan=False)


def test_credible_sets_and_scoping():
    left = pd.concat(
        [
            evidence(credible_set_id="cs1", pip=0.99),
            evidence(position=99, credible_set_id="cs1", pip=0.1),
        ]
    )
    right = pd.concat(
        [
            evidence("eqtl", credible_set_id="cs1", pip=0.97),
            evidence("eqtl", position=88, credible_set_id="cs1", pip=0.1),
            evidence("eqtl", gene="GENE2", credible_set_id="cs1", pip=0.5),
        ]
    )
    overlap = scrfu.tl.credible_set_overlap(left, right)
    assert len(overlap) == 2
    first = overlap.loc[overlap.right_set.str.contains('"GENE1"')].iloc[0]
    assert first.left_size == 2 and first.right_size == 2 and first.n_shared == 1
    assert first.n_shared_high_pip == 1 and first.jaccard == pytest.approx(1 / 3)
    assert json.loads(first.shared_variants) == ["1:12345:A:G"]
    result = scrfu.tl.regulatory_triangulation(left, eqtl=right)
    assert_frame_equal(result.credible_set_overlaps, overlap)
    assert result.matched_evidence.shared_high_pip.sum() == 1
    assert scrfu.tl.credible_set_overlap(evidence(), evidence("eqtl")).empty
    reverse = evidence("eqtl", ref="G", alt="A", effect_allele="A", credible_set_id="cs1")
    assert scrfu.tl.credible_set_overlap(left, reverse).empty


def test_sets_context_source_locus_and_missing_pip():
    left = evidence(credible_set_id="cs1")
    right = pd.concat(
        [
            evidence("eqtl", credible_set_id="cs1", context="CD4", locus_id="L1"),
            evidence("eqtl", credible_set_id="cs1", context="CD8", locus_id="L1"),
            evidence("eqtl", credible_set_id="cs1", context="CD4", locus_id="L2"),
        ]
    )
    overlap = scrfu.tl.credible_set_overlap(left, right)
    assert len(overlap) == 3 and overlap.n_shared_high_pip.eq(0).all()


def test_join_summary():
    result = scrfu.tl.regulatory_triangulation(evidence(), eqtl=evidence("eqtl"))
    summary = pd.DataFrame(
        dict(rfu_label=["RFU1", "unknown", "RFU1"], abundance=[3, 4, 5]), index=[9, 2, 8]
    )
    joined = scrfu.tl.join_regulatory_summary(summary, result)
    assert joined.index.tolist() == [9, 2, 8]
    assert joined.abundance.tolist() == [3, 4, 5]
    assert pd.isna(joined.loc[2, "has_eqtl"])
    with pytest.raises(ValueError, match="collision"):
        scrfu.tl.join_regulatory_summary(joined, result)
    with pytest.raises(ValueError, match="rfu_label"):
        scrfu.tl.join_regulatory_summary(pd.DataFrame(), result)


def test_plotting():
    mpl = pytest.importorskip("matplotlib")
    mpl.use("Agg")
    import matplotlib.pyplot as plt

    for raw in (evidence(), pd.DataFrame()):
        result = scrfu.tl.regulatory_triangulation(raw, eqtl=evidence("eqtl", pip=0.99))
        for plot in (scrfu.pl.regulatory_evidence_heatmap, scrfu.pl.regulatory_evidence_bar):
            fig, ax = plt.subplots()
            assert plot(result, ax=ax) is ax
            if plot is scrfu.pl.regulatory_evidence_heatmap and not raw.empty:
                np.testing.assert_array_equal(ax.images[0].get_array(), [[1, 0, 1, 0, 1]])
            plt.close(fig)


@pytest.fixture
def adapter():
    path = Path(__file__).parents[1] / "examples" / "matos_regulatory_triangulation.py"
    spec = importlib.util.spec_from_file_location("matos_example", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_matos_tensorqtl_and_ambiguity(adapter):
    raw = pd.DataFrame(
        [dict(variant_id="1:12345[b38]A,G", phenotype_id="GENE1", slope=0.4, pval_nominal=0.001)]
    )
    args = dict(
        layer="eqtl", format="tensorqtl", source="Matos", release="test", genome_build="hg38"
    )
    with pytest.raises(ValueError, match="allele_order"):
        adapter.adapt_matos_table(raw, **args)
    out = adapter.adapt_matos_table(raw, **args, allele_order="ref-alt")
    assert out.loc[0, "variant_key"] == "1:12345:A:G"
    assert out.loc[0, "gene"] == "GENE1" and out.loc[0, "pvalue"] == 0.001
    assert out.loc[0, "effect_allele"] is None
    with pytest.raises(ValueError, match="Embedded variant build"):
        adapter.adapt_matos_table(raw, **{**args, "genome_build": "hg19"}, allele_order="ref-alt")


def test_matos_susie_original_variant_and_coverage(adapter):
    raw = pd.DataFrame(
        [
            {
                "variant_id": "rs77",
                "variant_id.ss": "1:12345[b38]G,A",
                "peak": "peak1",
                "pip": 0.99,
                "cs": 1,
                "region": "1:10000-15000",
                "coverage": 0.10,
            }
        ]
    )
    out = adapter.adapt_matos_table(
        raw,
        layer="caqtl",
        format="susie",
        source="Matos",
        release="test",
        genome_build="hg38",
        allele_order="alt-ref",
    )
    assert out.loc[0, "variant_key"] == "1:12345:A:G"
    assert out.loc[0, "coverage"] == 0.10 and out.loc[0, "credible_set_id"] == "1"


def test_local_manifest_exports(adapter, tmp_path):
    evidence().to_csv(tmp_path / "rfu.tsv", sep="\t", index=False)
    pd.DataFrame([dict(variant_id="1:12345[b38]A,G", phenotype_id="GENE1", slope=0.2)]).to_csv(
        tmp_path / "eqtl.tsv", sep="\t", index=False
    )
    manifest = {
        "rfu_qtl": dict(
            source="user",
            release="v1",
            genome_build="hg38",
            files=[dict(path="rfu.tsv", format="canonical")],
        ),
        "eqtl": dict(
            source="Matos",
            release="fixture-v1",
            genome_build="hg38",
            context="CD4",
            files=[dict(path="eqtl.tsv", format="tensorqtl", allele_order="ref-alt")],
        ),
    }
    config = tmp_path / "manifest.json"
    config.write_text(json.dumps(manifest))
    adapter.run_manifest(config, tmp_path / "out")
    adapter.run_manifest(config, tmp_path / "out2")
    for path in (tmp_path / "out").iterdir():
        assert path.read_bytes() == (tmp_path / "out2" / path.name).read_bytes()
    prov = json.loads((tmp_path / "out" / "provenance.json").read_text())
    assert len(prov["inputs"]["eqtl"]["metadata"]["files"][0]["sha256"]) == 64
    hits = pd.read_csv(tmp_path / "out" / "regulatory_hits.tsv", sep="\t")
    assert len(hits) == 1 and hits.direction_concordant.isna().all()
    with pytest.raises(FileExistsError):
        adapter.run_manifest(config, tmp_path / "out")
    manifest["eqtl"]["files"][0]["sha256"] = "wrong"
    config.write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="SHA256"):
        adapter.run_manifest(config, tmp_path / "bad")


def test_position_precision_and_partial_coordinates():
    for value in ("9007199254740993", "12345.00000000000001"):
        with pytest.raises(ValueError, match="position"):
            scrfu.tl.regulatory_triangulation(evidence(position=value))
    raw = evidence(position=pd.NA, variant_id="1:12345:A:G")
    raw["position"] = raw.position.astype("Int64")
    assert scrfu.tl.regulatory_triangulation(raw).harmonized_rfu_qtl.loc[0, "position"] == 12345


def test_numeric_rfu_join_and_unmodified_index():
    result = scrfu.tl.regulatory_triangulation(evidence(rfu_label=1))
    summary = pd.DataFrame({"rfu_label": [1, 2, 1]}, index=[1, 1, 1])
    joined = scrfu.tl.join_regulatory_summary(summary, result)
    assert joined.rfu_label.tolist() == [1, 2, 1]
    assert joined.index.tolist() == [1, 1, 1]
    assert joined.has_rfu_qtl.iloc[0] and pd.isna(joined.has_rfu_qtl.iloc[1])


def test_multibase_variant_and_nonsignificant_presence():
    result = scrfu.tl.regulatory_triangulation(
        evidence(ref="AT", alt="A", effect_allele="A"),
        eqtl=evidence("eqtl", ref="AT", alt="A", effect_allele="A", pvalue=0.9),
    )
    assert result.matched_evidence.loc[0, "shared_variant_exact"]
    assert result.harmonized_rfu_qtl.loc[0, "has_eqtl"]


def test_gwas_invalid_effect_scales():
    for changes in ({"odds_ratio": 2}, {"beta": None, "odds_ratio": 0}):
        with pytest.raises(ValueError, match="beta or odds_ratio|odds_ratio must"):
            scrfu.tl.regulatory_triangulation(evidence(), gwas=evidence("gwas", **changes))


def test_matos_effect_column_and_canonical_map(adapter):
    raw = pd.DataFrame([dict(variant_id="1:12345[b38]A,G", phenotype_id="G", ea="G", slope=0.2)])
    out = adapter.adapt_matos_table(
        raw,
        layer="eqtl",
        format="tensorqtl",
        source="Matos",
        release="v1",
        genome_build="hg38",
        allele_order="ref-alt",
        effect_allele_column="ea",
    )
    assert scrfu.tl.regulatory_triangulation(evidence(), eqtl=out).matched_evidence.loc[
        0, "direction_concordant"
    ]
    renamed = adapter.adapt_matos_table(
        evidence("gwas").rename(columns={"beta": "effect"}),
        layer="gwas",
        format="canonical",
        source="GWAS",
        release="v1",
        genome_build="hg38",
        column_map={"effect": "beta"},
    )
    assert renamed.beta.iloc[0] == 0.4


def test_metadata_and_column_order_are_not_mutated():
    raw = evidence()
    meta = {"rfu_qtl": {"release": "v1", "source": "synthetic"}}
    before = json.dumps(meta)
    a = scrfu.tl.regulatory_triangulation(raw, input_metadata=meta)
    b = scrfu.tl.regulatory_triangulation(raw[raw.columns[::-1]], input_metadata=meta)
    assert a.harmonized_rfu_qtl.record_id.tolist() == b.harmonized_rfu_qtl.record_id.tolist()
    assert json.dumps(meta) == before


def test_regulatory_public_snapshot_is_experimental():
    snapshot = json.loads((Path(__file__).parents[1] / "docs" / "public_api_0.5.json").read_text())
    entries = [
        e
        for e in snapshot["entries"]
        if e["module"] == "scrfu.regulatory" or e["name"].startswith("regulatory_evidence_")
    ]
    assert len(entries) == 9
    assert all(e["stability"] == "experimental" for e in entries)
