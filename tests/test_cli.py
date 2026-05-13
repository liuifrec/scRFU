import pytest

from scrfu.cli import build_parser, main


def test_cli_parser():
    p = build_parser()
    ns = p.parse_args(
        [
            "call-rfu",
            "in.h5ad",
            "-o",
            "out.h5ad",
            "--rfu-dir",
            "/tmp/RFU",
            "--chain",
            "TRB",
            "--airr-key",
            "airr",
            "--out-key",
            "rfu",
            "--wrapper-r-path",
            "r/run_rfu_repo.R",
            "--rscript-bin",
            "/usr/bin/Rscript",
            "--extra-r-arg=--vanilla",
            "--extra-r-arg=--quiet",
        ]
    )
    assert ns.cmd == "call-rfu"
    assert ns.chain == "TRB"
    assert ns.rfu_dir == "/tmp/RFU"
    assert ns.wrapper_r_path == "r/run_rfu_repo.R"
    assert ns.rscript_bin == "/usr/bin/Rscript"
    assert ns.extra_r_arg == ["--vanilla", "--quiet"]


def test_cli_parser_requires_subcommand():
    parser = build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args([])


def test_cli_main_dispatches_current_call_rfu_signature(monkeypatch):
    calls: dict[str, object] = {}
    adata = object()

    def fake_read_h5ad(path: str) -> object:
        calls["input"] = path
        return adata

    def fake_call_rfu(adata_arg: object, **kwargs: object) -> None:
        calls["adata"] = adata_arg
        calls["kwargs"] = kwargs

    def fake_write_h5ad(adata_arg: object, path: str) -> None:
        calls["write"] = (adata_arg, path)

    monkeypatch.setattr("scrfu.cli.read_h5ad", fake_read_h5ad)
    monkeypatch.setattr("scrfu.cli.call_rfu", fake_call_rfu)
    monkeypatch.setattr("scrfu.cli.write_h5ad", fake_write_h5ad)

    main(
        [
            "call-rfu",
            "input.h5ad",
            "-o",
            "output.h5ad",
            "--rfu-dir",
            "/rfu",
            "--wrapper-r-path",
            "r/run_rfu_repo.R",
            "--rscript-bin",
            "Rscript",
            "--workdir",
            ".scrfu/test-run",
            "--extra-r-arg=--vanilla",
        ]
    )

    assert calls["input"] == "input.h5ad"
    assert calls["adata"] is adata
    assert calls["write"] == (adata, "output.h5ad")
    assert calls["kwargs"] == {
        "rfu_dir": "/rfu",
        "chain": "TRB",
        "airr_key": "airr",
        "out_key": "rfu",
        "wrapper_r_path": "r/run_rfu_repo.R",
        "rscript_bin": "Rscript",
        "extra_r_args": ["--vanilla"],
        "workdir": ".scrfu/test-run",
    }
