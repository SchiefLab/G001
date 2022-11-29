from click.testing import CliRunner
from g001.app import main
from pytest import TempdirFactory


def test_process(tmpdir_factory: TempdirFactory) -> None:
    """test collation directive"""
    tmp_output = tmpdir_factory.mktemp("test_process")
    runner = CliRunner()
    for x in ["vrc", "fhcrc"]:
        result = runner.invoke(
            main,
            [
                "process-flow",
                "--verbose",
                "-s",
                x,
                "-i",
                "flow_input_sampled/",
                "-o",
                str(tmp_output),
            ],
        )
        assert result.exit_code == 0
