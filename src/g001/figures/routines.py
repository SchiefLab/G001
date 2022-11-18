# -*- coding: utf-8 -*-
import json
import logging
import re
from pathlib import Path
from typing import List, Union

import pandas as pd
from sadie.airr import LinkedAirrTable

logger = logging.getLogger(__name__)


def get_iter_kabat(x):
    """only get VH/VK/VL region"""
    return_list = []
    if isinstance(x, str):
        x = eval(x)
    for y in x:
        number = re.findall(r"\d+", y)[0]
        if int(number) > 93:
            continue
        return_list.append(y)
    return return_list


def tag_vrc01_class(dataframe: LinkedAirrTable, v_gene: str = "IGHV1-2*02") -> LinkedAirrTable:
    """Given a LinkedAirrTable, tag all vrc01 class antibodies. Might have to change v_gene if you are using
    chimerized models and the v_gene is tagged 'human|IGHV1-2*02'

    Parameters
    ----------
    dataframe : LinkedAirrTable
        The input DataFrame
    v_gene : str, optional
        the heavy v_gene to look for, by default "IGHV1-2*02"

    Returns
    -------
    LinkedAirrTable
        A linked airrtable with is_vrc01_class field added
    """
    if not isinstance(dataframe, LinkedAirrTable):
        raise TypeError(f"{type(dataframe)} must be LinkedAirrTable")
    vrc01_class_index = dataframe[
        (dataframe["v_call_heavy"].str.contains(v_gene)) & (dataframe["cdr3_aa_light"].str.len() == 5)
    ].index

    if vrc01_class_index.empty:
        logger.warning(f"Warning, no VRC01 class found, maybe {v_gene} is not the right call")
    dataframe.loc[:, "is_vrc01_class"] = False
    dataframe.loc[vrc01_class_index, "is_vrc01_class"] = True
    return LinkedAirrTable(dataframe)


def frequency_dataframe(dataframe: pd.DataFrame, groupby: Union[str, List[str]]) -> pd.DataFrame:
    """Generate a frequency dataframe based on a groupby condition

    Parameters
    ----------
    dataframe : pd.DataFrame
        The intiial dataframe to make a frequency_dataframe
    groupby : Union[str, List[str]]
        The field or list of fields to group on


    Returns
    -------
    dataframe
        The frequency dataframe
    """
    if "is_vrc01_class" not in dataframe.columns:
        raise ValueError("is_vrc01_class must be a field in input dataframe")
    grouped_df = []
    groupby_immunization = dataframe.groupby(groupby)
    for group, group_df in groupby_immunization:
        new_entry = {
            "vrc01_class": len(group_df.query("is_vrc01_class==True")),
            "total_pairs": len(group_df),
        }
        new_entry["vrc01_frequency"] = (new_entry["vrc01_class"] / new_entry["total_pairs"]) * 100
        if isinstance(groupby, list):
            for x, y in zip(group, groupby):
                new_entry[y] = x
            columns_list = groupby + ["vrc01_class", "total_pairs", "vrc01_frequency"]
        else:
            new_entry[groupby] = group
            columns_list = [groupby] + ["vrc01_class", "total_pairs", "vrc01_frequency"]
        grouped_df.append(pd.Series(new_entry))

    # frequency df will be our plottable dataframe
    frequency_df = pd.DataFrame(grouped_df)[columns_list]
    return frequency_df


def find_100b(row):
    "If we have a tryptophan at the purported 100b position, then return true"
    if not isinstance(row, str) or not row:
        return False
    elif len(row) < 6:
        return False
    elif row[-6] == "W":
        return True
    return False


def add_mutational_sets(
    dataframe: LinkedAirrTable,
    vh12_reference_airr_table: pd.DataFrame,
    cottrell_path: Path
):
    if not isinstance(dataframe, LinkedAirrTable):
        try:
            dataframe = LinkedAirrTable(dataframe)
        except Exception:
            raise TypeError(f"{dataframe} must be type of LinkedAirrTable or castable to linkedAirrTable")
    if "is_vrc01_class" not in dataframe.columns:
        raise ValueError(f"is_vrc01_class field not found in {dataframe.columns}")
    if "mutations_heavy" not in dataframe.columns:
        raise ValueError(
            f"""mutations_heavy field not found in {dataframe.columns}.
            Run sadie.airr methods `run_mutational_analysis` to generate"""
        )
    # vh12_reference_airr_table = pd.read_feather(
    #     Path(__file__).parent.joinpath("data/vrc1-2_mabs_extended_airr.feather")
    # )
    cottrell_super_focus = json.load(open(cottrell_path))
    cottrell_super_focus_positive = cottrell_super_focus["positive_set"]
    cottrell_super_focus_negative = cottrell_super_focus["negative_set"]

    cottrell_mabs = "N6 VRC27 VRC01 12A12 PCIN63_71I VRC-PG20".split()
    jardine_mabs = "12A12 3BNC60 VRC-PG04 VRC-PG20 VRC-CH31 VRC01".split()
    cotrell_focus = vh12_reference_airr_table[vh12_reference_airr_table["sequence_id"].isin(cottrell_mabs)]
    jardine_focus = vh12_reference_airr_table[vh12_reference_airr_table["sequence_id"].isin(jardine_mabs)]

    cotrell_focus_heavy_sets = set([item for sublist in cotrell_focus["mutations_heavy"].to_list() for item in sublist])
    cotrell_focus_light_sets = set([item for sublist in cotrell_focus["mutations_light"].to_list() for item in sublist])

    jardine_focus_heavy_sets = set([item for sublist in jardine_focus["mutations_heavy"].to_list() for item in sublist])
    jardine_focus_light_sets = set([item for sublist in jardine_focus["mutations_light"].to_list() for item in sublist])

    dataframe["cottrell_focused_v_common_heavy_positive"] = dataframe["mutations_heavy"].apply(
        lambda x: list(set(get_iter_kabat(x)).intersection(cottrell_super_focus_positive))
    )
    dataframe["cottrell_focused_v_common_heavy_negative"] = dataframe["mutations_heavy"].apply(
        lambda x: list(set(map(lambda y: y[1:-1], x)).intersection(cottrell_super_focus_negative))
    )
    dataframe["cottrell_focused_v_common_score"] = dataframe["cottrell_focused_v_common_heavy_positive"].apply(
        lambda x: len(x)
    ) - dataframe["cottrell_focused_v_common_heavy_negative"].apply(lambda x: len(x))

    dataframe["100bW"] = dataframe["junction_aa_heavy"].apply(find_100b)
    dataframe["cottrell_focused_v_common_score"] += dataframe["100bW"].apply(lambda x: {True: 1, False: 0}[x])

    dataframe["cottrell_v_common_heavy"] = dataframe["mutations_heavy"].apply(
        lambda x: list(set(get_iter_kabat(x)).intersection(cotrell_focus_heavy_sets))
    )
    dataframe["cottrell_v_common_heavy_score"] = dataframe["cottrell_v_common_heavy"].apply(lambda x: len(x))

    dataframe["cottrell_v_common_light"] = dataframe["mutations_light"].apply(
        lambda x: list(set(get_iter_kabat(x)).intersection(cotrell_focus_light_sets))
    )
    dataframe["cottrell_v_common_light_score"] = dataframe["cottrell_v_common_light"].apply(lambda x: len(x))

    dataframe["jardine_v_common_heavy"] = dataframe["mutations_heavy"].apply(
        lambda x: list(set(get_iter_kabat(x)).intersection(jardine_focus_heavy_sets))
    )
    dataframe["jardine_v_common_heavy_score"] = dataframe["jardine_v_common_heavy"].apply(lambda x: len(x))
    dataframe["jardine_v_common_light"] = dataframe["mutations_light"].apply(
        lambda x: list(set(get_iter_kabat(x)).intersection(jardine_focus_light_sets))
    )
    dataframe["jardine_v_common_light_score"] = dataframe["jardine_v_common_light"].apply(lambda x: len(x))
    return dataframe
