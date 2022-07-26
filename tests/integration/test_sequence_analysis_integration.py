from functools import cache
from click.testing import CliRunner
from g001.app import click_run_sa
from g001.data import Data
from pytest import TempdirFactory, fixture
from pathlib import Path
import pandas as pd
from glob import glob
from g001.sequence_pipeline import (
    find,
    model,
    ngs,
    split,
    correct,
    annotate,
    join,
    tag,
    swap,
    pair,
    unblind,
    personalize,
    mutation,
)
import pytest


@cache
@pytest.fixture
def tmp_path(tmp_path_factory: pytest.TempPathFactory) -> Path:
    tmp_path = tmp_path_factory.mktemp("test_sequence_analysis_integration")
    # return Path("tests/output/")
    return Path(tmp_path)


@pytest.fixture(autouse=True)
def data() -> Data:
    data = Data("./data")
    return data


@cache
def grab_fixture_df(name: str) -> pd.DataFrame:
    if name in ["pair", "unblind", "personalize", "mutation"]:
        df = pd.read_feather(f"tests/fixtures/{name}.feather")
        return df.astype(
            {
                "d_alignment_start_light": float,
                "d_alignment_end_light": float,
                "d_sequence_start_light": float,
                "d_sequence_end_light": float,
                "np2_length_light": float,
                "d_germline_start_light": float,
                "d_germline_end_light": float,
            }
        )
    else:
        return pd.read_feather(f"tests/fixtures/{name}.feather")


def assert_frames_equal(tmp_path: Path, current: pd.DataFrame, fixture_df: pd.DataFrame, test_name: str) -> None:
    out_path = tmp_path / f"{test_name}.feather"
    print("Writing to:", out_path)
    current.to_feather(out_path)
    print("Reading", out_path)
    current = pd.read_feather(out_path)
    pd.testing.assert_frame_equal(current, fixture_df)


def test_find_module(tmp_path: Path, data: Data) -> None:
    working_dataframe = find(data)
    assert_frames_equal(tmp_path, working_dataframe, grab_fixture_df("find"), "find")


def test_model_module(tmp_path: Path, data: Data) -> None:
    working_dataframe = grab_fixture_df("find")
    working_dataframe = model(data, working_dataframe)
    assert_frames_equal(tmp_path, working_dataframe, grab_fixture_df("model"), "model")


def test_split_module(tmp_path: Path, data: Data) -> None:
    working_dataframe = grab_fixture_df("model")
    working_dataframe = split(data, working_dataframe)
    assert_frames_equal(tmp_path, working_dataframe, grab_fixture_df("split"), "split")


def test_correct_module(tmp_path: Path, data: Data) -> None:
    working_dataframe = grab_fixture_df("split")
    working_dataframe = correct(data, working_dataframe)
    assert_frames_equal(tmp_path, working_dataframe, grab_fixture_df("correct"), "correct")


@pytest.mark.skip(reason="Not implemented yet")
def test_annotate_module(tmp_path: Path, data: Data) -> None:
    working_dataframe = grab_fixture_df("correct")
    working_dataframe = annotate(data, working_dataframe)
    assert_frames_equal(tmp_path, working_dataframe, grab_fixture_df("annotate"), "annotate")


def test_join_module(tmp_path: Path, data: Data) -> None:
    working_dataframe = grab_fixture_df("annotate")
    working_dataframe = join(data, working_dataframe)
    assert_frames_equal(tmp_path, working_dataframe, grab_fixture_df("join"), "join")


def test_ngs_module(tmp_path: Path, data: Data) -> None:
    working_dataframe = grab_fixture_df("join")
    working_dataframe = ngs(data, working_dataframe)
    assert_frames_equal(tmp_path, working_dataframe, grab_fixture_df("ngs"), "ngs")


def test_tag_module(tmp_path: Path, data: Data) -> None:
    working_dataframe = grab_fixture_df("ngs")
    working_dataframe = tag(data, working_dataframe)
    assert_frames_equal(tmp_path, working_dataframe, grab_fixture_df("tag"), "tag")


def test_swap_module(tmp_path: Path, data: Data) -> None:
    working_dataframe = grab_fixture_df("tag")
    working_dataframe = swap(data, working_dataframe)
    assert_frames_equal(tmp_path, working_dataframe, grab_fixture_df("swap"), "swap")


def test_pair_module(tmp_path: Path, data: Data) -> None:
    working_dataframe = grab_fixture_df("swap")
    working_dataframe, candidates_df = pair(data, working_dataframe)
    assert_frames_equal(tmp_path, working_dataframe, grab_fixture_df("pair"), "pair")
    assert_frames_equal(tmp_path, candidates_df, grab_fixture_df("pairing_candidates"), "pairing_candidates")


def test_unblind_module(tmp_path: Path, data: Data) -> None:
    working_dataframe = grab_fixture_df("pair")
    working_dataframe = unblind(data, working_dataframe)
    assert_frames_equal(tmp_path, working_dataframe, grab_fixture_df("unblind"), "unblind")


@pytest.mark.skip(reason="Not implemented yet")
def test_personalize_module(tmp_path: Path, data: Data) -> None:
    working_dataframe = grab_fixture_df("unblind")
    working_dataframe = personalize(data, working_dataframe)
    assert_frames_equal(tmp_path, working_dataframe, grab_fixture_df("personalize"), "personalize")


@pytest.mark.skip(reason="Not implemented yet")
def test_mutate_module(tmp_path: Path, data: Data) -> None:
    working_dataframe = grab_fixture_df("personalize")
    working_dataframe = mutation(data, working_dataframe)
    assert_frames_equal(tmp_path, working_dataframe, grab_fixture_df("mutation"), "mutation")
