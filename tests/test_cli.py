import pandas as pd
from scrfu.cli import build_parser


def test_cli_parser():
    p = build_parser()
    ns = p.parse_args(
        [
            "call-rfu",
            "in.h5ad",
            "-o",
            "out.h5ad",
            "--chain",
            "TRB",
            "--airr-key",
            "airr",
            "--out-key",
            "rfu",
            "--rscript",
            "r/RFU_call.R",
            "--atlas",
            "data/km5000_centroids.tsv.gz",
        ]
    )
    assert ns.cmd == "call-rfu"
    assert ns.chain == "TRB"