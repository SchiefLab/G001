from click.testing import CliRunner
from g001.app import main


def test_collate():
    """test collation directive"""
    runner = CliRunner()
    result = runner.invoke(
        main,
        [
            "collate",
            "-1",
            "tests/data/manifests/FHCRC_Flow_Manifest.csv",
            "-2",
            "tests/data/manifests/VRC_Flow_Manifest.csv",
            "-f",
            "data/flow/flow_processed_out",
            "-o",
            "collated_flow"
        ],
    )
    assert result.exit_code == 0
