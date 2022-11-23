from click.testing import CliRunner
from g001.app import main
from pytest import TempdirFactory


def test_combine(tmp_path_factory: TempdirFactory) -> None:
    """test collation directive"""
    path = tmp_path_factory.mktemp("test_combine")
    print(path)
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
            "data/flow/collated_flow/",
            "-o",
            "combined_data",
        ],
    )
    if result.exit_code != 0:
        print(result.output)
        print(result.stderr)
    print(result.stdout)
    assert result.exit_code == 0
