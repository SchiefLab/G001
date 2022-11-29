from click.testing import CliRunner
from g001.app import main
from pytest import TempdirFactory


def test_combine(tmp_path_factory: TempdirFactory) -> None:
    """test collation directive"""
    path = tmp_path_factory.mktemp("test_combine")
    print(path)
    runner = CliRunner()
    args = ["combine", "--verbose", "-o", str(path)]
    print("Running", "g001 combine", " ".join(args))
    result = runner.invoke(main, args)
    print(result.stdout)
    assert result.exit_code == 0
