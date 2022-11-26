from click.testing import CliRunner
from g001.app import main
from pytest import TempdirFactory


def test_collate(tmpdir_factory: TempdirFactory) -> None:
    """test collation directive"""
    tmp_output = tmpdir_factory.mktemp("test_collate")
    runner = CliRunner()
    result = runner.invoke(
        main,
        [
            "collate",
            "--verbose",
            "-o",
            str(tmp_output),
        ],
    )
    assert result.exit_code == 0
    result = runner.invoke(
        main,
        [
            "collate",
            "--verbose",
            "-1",
            "data/flow/flow_processed_out/fhcrc/fhcrc_manifest.csv",
            "-2",
            "data/flow/flow_processed_out/vrc/vrc_manifest.csv",
            "-f",
            "data/flow/flow_processed_out",
            "-s",
            "data/flow/swap.csv",
            "-o",
            str(tmp_output),
        ],
    )
    assert result.exit_code == 0
