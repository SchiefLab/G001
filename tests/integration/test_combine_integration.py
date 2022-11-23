from click.testing import CliRunner
from g001.app import main


def test_collate():
    """test collation directive"""
    runner = CliRunner()
    result = runner.invoke(
        main, ["combine", "--verbose", "-s", "data/sequence/", "-c", "data/flow/collated_flow/", "-o", "combined_data/"]
    )
    if result.exit_code != 0:
        print(result.output)
        print(result.stderr)
    assert result.exit_code == 0
