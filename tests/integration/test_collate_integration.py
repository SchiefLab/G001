from click.testing import CliRunner
from g001.app import main


def test_collate():
    """test collation directive"""
    runner = CliRunner()
    result = runner.invoke(main, ["collate", "-f", "data/flow/flow_processed_out"])
    assert result.exit_code == 0
