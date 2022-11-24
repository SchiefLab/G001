from click.testing import CliRunner
from g001.app import main
from pytest import TempdirFactory


def test_combine(tmp_path_factory: TempdirFactory) -> None:
    """test collation directive"""
    path = tmp_path_factory.mktemp("test_combine")
    print(path)
    runner = CliRunner()
    args = [
        "combine",
        "--verbose",
        "-1",
        "data/flow/flow_processed_out/fhrc/fhcrc_manifest.csv",
        "-2",
        "data/flow/flow_processed_out/vrc/vrc_manifest.csv",
        "-s",
        "data/sequence/",
        "-c",
        "data/flow/collated_flow/",
        "-o",
        "combined_data",
    ]
    print("Running", "g001 combine", " ".join(args))
    result = runner.invoke(main, args)
    print(result.stdout)
    assert result.exit_code == 0
