from click.testing import CliRunner
import pytest
from g001.app import main


@pytest.mark.order(2)
def test_combine():
    """test collation directive"""
    runner = CliRunner()
    result = runner.invoke(
        main,
        [
            "combine",
            "--verbose",
            "-1",
            "tests/data/manifests/FHCRC_Flow_Manifest.csv",
            "-2",
            "tests/data/manifests/VRC_Flow_Manifest.csv",
            "-s",
            "data/sequence/",
            "-c",
            "data/flow/collated_flow/",  # TODO: this will eventually point to g001 collate output
            "-o",
            "tests/data/flow/combined_data/",
        ],
    )
    if result.exit_code != 0:
        print(result.output)
        print(result.stderr)
    assert result.exit_code == 0
