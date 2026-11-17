#!/usr/bin/env python3
"""Synthetic Scirpy/MuData scRFU tutorial; mock assignments are not scientific results."""

from __future__ import annotations

import argparse
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
import scanpy as sc
import scirpy as ir
from anndata import AnnData
from mudata import MuData, read_h5mu
from scirpy.io import AirrCell, from_airr_cells

import scrfu
from scrfu.rfu import RFURunResult


def _synthetic_run(self: object, features: pd.DataFrame, **kwargs: object) -> RFURunResult:
    del self, kwargs
    rows = features.copy().reset_index(drop=True)
    codes, _ = pd.factorize(rows["cdr3aa"], sort=True)
    rows["eligibility_status"] = "eligible"
    rows["unique_sequence_id"] = [f"synthetic_sequence_{code:04d}" for code in codes]
    rows["rfu_id"] = pd.array(codes + 1, dtype="Int64")
    rows["rfu_label"] = pd.Series([f"mock_RFU{code + 1}" for code in codes], dtype="string")
    rows["rfu_score"] = 0.9
    rows["pass_thr"] = pd.array([True] * len(rows), dtype="boolean")
    rows["rfu_status"] = "assigned_threshold_pass"
    return RFURunResult(
        rows,
        "",
        "",
        0,
        metadata={"rfu_threshold": 0.6, "synthetic_mock": True},
    )


def _external_files(root: Path) -> tuple[Path, Path]:
    rfu_dir = root / "synthetic-not-official-rfu"
    rfu_dir.mkdir(parents=True, exist_ok=True)
    (rfu_dir / "RFU.R").write_text("AssignRFUs <- function(ff) {}\n", encoding="utf-8")
    (rfu_dir / "trimerMDSfit_small.Rdata").write_bytes(b"synthetic trimer")
    (rfu_dir / "km5000noMax.Rdata").write_bytes(b"synthetic centroids")
    wrapper = root / "synthetic_wrapper.R"
    wrapper.write_text("# synthetic tutorial placeholder\n", encoding="utf-8")
    return rfu_dir, wrapper


def _airr() -> AnnData:
    cells = []
    for index, cdr3 in enumerate(("CASSLG", "CASSLG", "CASSQG", "CASSPG")):
        cell = AirrCell(f"cell-{index}")
        for locus, junction, v_call, j_call in (
            ("TRA", f"CAV{index}F", "TRAV1", "TRAJ1"),
            ("TRB", cdr3, f"TRBV{index % 2 + 1}", "TRBJ1-1"),
        ):
            chain = AirrCell.empty_chain_dict()
            chain.update(
                {
                    "locus": locus,
                    "junction": "TGTGCC",
                    "junction_aa": junction,
                    "v_call": v_call,
                    "j_call": j_call,
                    "productive": True,
                }
            )
            cell.add_chain(chain)
        cells.append(cell)
    adata = from_airr_cells(cells)
    adata.obs["sample"] = ["s1", "s1", "s2", "s2"]
    adata.obs["cell_type"] = ["CD8 T", "CD8 T", "CD4 T", "CD4 T"]
    adata.obs["cell_state"] = ["memory", "effector", "naive", "memory"]
    adata.obsm["X_umap"] = np.asarray([[0, 0], [1, 0], [0, 1], [1, 1]], dtype=float)
    return adata


def run(outdir: Path) -> Path:
    outdir.mkdir(parents=True, exist_ok=True)
    airr = _airr()
    ir.pp.index_chains(airr)
    gex = AnnData(X=np.eye(airr.n_obs), obs=airr.obs.copy())
    gex.obsm["X_umap"] = airr.obsm["X_umap"].copy()
    gex.obsm["protein"] = np.arange(airr.n_obs * 2).reshape(airr.n_obs, 2)
    mdata = MuData({"airr": airr, "gex": gex})

    rfu_dir, wrapper = _external_files(outdir)
    with patch("scrfu.backends.rfu_repo.RFURepoBackend.run", _synthetic_run):
        scrfu.tl.assign_rfu(
            mdata,
            rfu_dir=rfu_dir,
            wrapper_r_path=wrapper,
            cell_summary=True,
            summary_policy="ambiguity_aware",
        )

    assert scrfu.tl.validate_scrfu_schema(mdata)["status"] == "ok"
    assigned = mdata.mod["airr"].obs.rename(
        columns={
            "scrfu_rfu_label": "rfu_label",
            "scrfu_threshold_pass": "pass_thr",
        }
    )
    assigned["cell_id"] = assigned.index.astype(str)
    assigned["chain"] = "TRB"
    assigned["cdr3aa"] = ["CASSLG", "CASSLG", "CASSQG", "CASSPG"]
    scrfu.tl.rfu_pseudobulk(assigned, sample_key="sample")
    scrfu.tl.rfu_phenotype_coupling(assigned, phenotype_key="cell_type", sample_key="sample")
    sc.pl.umap(mdata.mod["airr"], color="scrfu_rfu_label", show=False)

    output = outdir / "synthetic_scverse_tutorial.h5mu"
    mdata.write_h5mu(output)
    loaded = read_h5mu(output)
    scrfu.tl.validate_scrfu_schema(loaded)
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    output = run(args.outdir.resolve())
    print(f"Wrote {output}")
    print("All RFU values in this tutorial are synthetic mock assignments.")


if __name__ == "__main__":
    main()
