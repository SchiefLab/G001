# -*- coding: utf-8 -*-
from typing import Any, Dict, List, Optional, Tuple
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt


def get_stats(df, kd_lookup="KD_fix"):
    """Get table stats, counts, median and donors represented"""
    _count = int(len(df))
    less_than_100uM = (len(df[df[kd_lookup] < 1e-4]) / _count) * 100
    median = df[kd_lookup].median()
    subjects = int(len(df["PubID"].unique()))
    return pd.Series({"count": _count, "greater_than_100": less_than_100uM, "median": median, "num_subjects": subjects})


def plot_binding_table(
    ax: plt.Axes,
    spr_for_plot: pd.DataFrame,
    annotate_args: Dict[Any, Any],
    title: str,
    title_xy: Tuple[float, float],
    plot_row_label: bool = False,
    overide_stats: Optional[int] = None,
    groupby_col: List[str] = ["weeks_post", "is_vrc01_class"],
    fontsize=8,
    table_scale=(0.97, 1),
    hue_order=None,
    median_label="Median $\mathregular{K_D}$ (nM)",  # noqa: W605
    median_scale=1e9,
    ascending_rank=[1, 0],
):
    if plot_row_label:
        row_label = [
            "N Abs tested",
            "% $\mathregular{K_D}$s < 100 $\mu$M",  # noqa: W605
            median_label,
            "N Participants",
        ]
    else:
        row_label = None
    stats = (
        spr_for_plot.groupby(groupby_col)
        .apply(lambda x: get_stats(x))
        .transpose()
        .sort_index(level=[0, 1], axis=1, ascending=ascending_rank)
    )

    if overide_stats:
        stats.loc["num_subjects", :] = overide_stats

    # convert
    gt_median = r"$\mathregular{-}$"
    # gt_median = r"$\mathregular{\geq100}$"
    cell_text = []
    nabs = list(stats.values[0])
    gt_100 = list(stats.values[1])
    median = list(stats.values[2])
    donors = list(stats.values[3])

    def scale_median(x):
        if np.isnan(x):
            return "-"
        if x != 1e-4:
            i = median_scale * x
            if i < 1:
                return round(i, 1)
            else:
                return int(i)
        else:
            return gt_median

    nabs = [int(i) if not np.isnan(i) else "-" for i in nabs]
    gt_100 = [int(i) if not np.isnan(i) else "-" for i in gt_100]
    median = [scale_median(i) for i in median]
    donors = [int(i) if not np.isnan(i) else "-" for i in donors]
    cell_text = [nabs, gt_100, median, donors]
    if hue_order:
        new_cell_text = []
        for row in cell_text:
            i = len(hue_order)
            while i < len(row):
                row.insert(i, "")
                i += len(hue_order) + 1
            new_cell_text.append(row)
        cell_text = new_cell_text
    table = ax.table(
        cellText=cell_text,
        cellLoc="center",
        rowLoc="right",
        rowLabels=row_label,
        edges="open",
        loc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(fontsize)
    table.scale(table_scale[0], table_scale[1])
    ax.annotate(title, xy=title_xy, **annotate_args)
    ax.axis("off")


def plot_core_table(
    axis_apec: plt.axis,
    axis_df: plt.axis,
    annotate_args: dict,
    plot_row_titles: bool,
    title: str,
    hue_order: list,
    kd_lookup: str = "KD_fix",
    override_subjects: bool = False,
):
    if plot_row_titles:
        row_label = [
            "N Abs tested",
            "% $\mathregular{K_D}$s < 100 $\mu$M",  # noqa: W605
            "Median $\mathregular{K_D}$ (nM)",  # noqa: W605
            "N Participants",
        ]
    else:
        row_label = None

    axis_df["Analyte"] = pd.Categorical(axis_df["Analyte"], categories=hue_order)
    stats = (
        axis_df.sort_values("Analyte").groupby(["x-axis", "Analyte"]).apply(lambda x: get_stats(x, kd_lookup=kd_lookup))
    )
    if override_subjects:
        # if CLK/HuGL, input 12
        stats["num_subjects"] = 12
    stats = stats.transpose()
    stats = stats[sorted(stats.columns, key=lambda x: int(x[0].split("\n")[0].split()[1]))]
    # get all cell stats and get them back as a List[List] of values

    def intergize(y):
        if y == "-" or y == "gt":
            return r"$\mathregular{\geq 10^5}$"
        try:
            return round(y)
        except ValueError:
            return y

    cell_text = list(
        map(
            lambda x: list(map(lambda y: intergize(y), x)),
            stats.fillna(0).values,
        )
    )
    new_cell_text = []
    for row in cell_text:
        i = len(hue_order)
        while i < len(row):
            row.insert(i, "")
            i += len(hue_order) + 1
        new_cell_text.append(row)
    table = axis_apec.table(
        cellText=cell_text,
        cellLoc="center",
        rowLoc="right",
        rowLabels=row_label,
        edges="open",
        loc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(0.97, 1)
    xy = (0.5, 1.8)
    axis_apec.annotate(title, xy=xy, **annotate_args)
    axis_apec.axis("off")


def plot_boost_table(
    axis_apec: plt.axis,
    axis_df: plt.axis,
    annotate_args: dict,
    plot_row_titles: bool,
    title: str,
    hue_order: list,
    kd_lookup: str = "KD_fix",
    override_subjects: bool = False,
):
    if plot_row_titles:
        row_label = [
            "N Abs tested",
            "% $\mathregular{K_D}$s < 100 $\mu$M",  # noqa: W605
            "Median $\mathregular{K_D}$ (nM)",  # noqa: W605
            "N Participants",
        ]
    else:
        row_label = None

    axis_df["PublishAnalyte"] = pd.Categorical(axis_df["PublishAnalyte"], categories=hue_order)
    stats = (
        axis_df.sort_values("PublishAnalyte")
        .groupby(["x-axis", "PublishAnalyte"])
        .apply(lambda x: get_stats(x, kd_lookup=kd_lookup))
    )
    if override_subjects:
        # if CLK/HuGL, input 12
        stats["num_subjects"] = 12
    stats = stats.transpose()
    stats = stats[sorted(stats.columns, key=lambda x: int(x[0].split("\n")[0].split()[1]))]
    # get all cell stats and get them back as a List[List] of values

    def intergize(y):
        if y == "-" or y == "gt":
            return r"$\mathregular{\geq 10^5}$"
        try:
            return round(y)
        except ValueError:
            return y

    cell_text = list(
        map(
            lambda x: list(map(lambda y: intergize(y), x)),
            stats.fillna(0).values,
        )
    )
    new_cell_text = []
    for row in cell_text:
        i = len(hue_order)
        while i < len(row):
            row.insert(i, "")
            i += len(hue_order) + 1
        new_cell_text.append(row)
    table = axis_apec.table(
        cellText=cell_text,
        cellLoc="center",
        rowLoc="right",
        rowLabels=row_label,
        edges="open",
        loc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(0.97, 1)
    xy = (0.5, 2.0)
    axis_apec.annotate(title, xy=xy, **annotate_args)
    axis_apec.axis("off")
