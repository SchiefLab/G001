from __future__ import annotations
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

from g001.data import Data
from g001.figures.util import adjust_boxplot, adjusts_axis, MinorSymLogLocator
import matplotlib.ticker as mtick
import math


def plot_response_panel(
    row_axis: plt.axes,
    plottable_df: pd.DataFrame,
    combined_timepoints: dict,
    palette: dict,
    annotate_args: dict,
    y_min: int,
    y_max: int,
    x_value: str = "Treatment",
    y_value: str = "Response",
    y_label: str = "% VRC01-class response\namong participant",
    labelpad: int = 4,
    plot_x_labels: bool = True,
    skip_header: bool = True,
):
    """
    Logic that handles resonse counts per timepoint - usually panel C
    """
    # top row panel A is VRC01 class counts
    timepoints = list(combined_timepoints.values())
    if isinstance(y_label, str):
        y_label = y_label
    elif isinstance(y_label, list) and len(y_label) == 3:
        y_label_pbmc, y_label_gc, y_label_pb = y_label
        y_label_pad_pbmc, y_label_pad_gc, y_label_pad_pb = labelpad
    else:
        raise ValueError("y_label must be a string or a list of 3 strings")

    # keep a counter since it won't enumerate when we skip the dummies
    weeks_index = 0
    for ax_index, ax in enumerate(row_axis):
        # remove the dummy axis
        if ax_index == 5 or ax_index == 8:
            ax.remove()
            continue

        # week lookup is interger of week
        week = timepoints[weeks_index]

        # Plotting df is only the week we are interested in
        plotting_df = plottable_df[plottable_df["weeks_post"] == week]
        sns.pointplot(ax=ax, data=plotting_df, x=x_value, y=y_value, palette=palette, ci="wilson")

        weeks_index += 1
        ax.set_ylim(y_min, y_max)
        ax = adjusts_axis(
            ax, ax_index, week, y_label, annotate_args, plot_x_labels=plot_x_labels, skip_header=skip_header
        )

        # add the ticks at every 0.1
        ax.yaxis.set_major_locator(mtick.MultipleLocator(base=0.1))

        # take the 0.1 and X by 100
        ax.yaxis.set_major_formatter(mtick.FuncFormatter(lambda y, _: str(int(float(y) * 100))))
        if ax_index > 0:
            ax.set_yticklabels([])
        if ax_index == 0 and isinstance(y_label, str):
            ax.set_ylabel(y_label, labelpad=labelpad)
        elif ax_index in [0, 6, 9] and isinstance(y_label, list):
            if ax_index == 0:
                ax.set_ylabel(y_label_pbmc, labelpad=y_label_pad_pbmc)
            elif ax_index == 6:
                ax.set_ylabel(y_label_gc, labelpad=y_label_pad_gc)
            else:
                ax.set_ylabel(y_label_pb, labelpad=y_label_pad_pb)
    return row_axis


def plot_frequency_panel(
    row_axis: plt.axes,
    plottable_df: pd.DataFrame,
    combined_timepoints: dict,
    palette: dict,
    annotate_args: dict,
    y_min: int,
    y_max: int,
    x_value: str = "Treatment",
    y_value: str | list[str] = "Percent of IgG+ B cells detected as VRC01-class",
    y_label: str | list[str] = "% VRC01-class\namong IgG+ B cells",
    thresh: float = 0.0001,
    baseline_jitter: int = 0.35,
    baseline_alpha: int = 0.9,
    value_jitter: int = 0.1,
    value_alpha: int = 0.9,
    strip_size: int = 4,
    jitter_buff: int = 0.15,
    skip_header: bool = True,
    y_text_pos: int | None = 20,
    x_tick_labels: bool = False,
    scale: str = "symlog",
    labelpad: int | list[int] = 0,
    plot_placebos_seperate: bool = True,
    replace_0_with_nan: bool = False,
) -> plt.figure:
    """
    Logic that handles frequency per timepoint - usually panel B
    """
    # top row panel A is VRC01 class counts
    timepoints = list(combined_timepoints.values())

    if isinstance(y_label, str):
        y_label = y_label
    elif isinstance(y_label, list) and len(y_label) == 3:
        if not isinstance(labelpad, list):
            raise Exception("if you pass a list of y labels, you need a list of labelpads, eg [0, 0, 0]")
        y_label_pbmc, y_label_gc, y_label_pb = y_label
        y_label_pad_pbmc, y_label_pad_gc, y_label_pad_pb = labelpad
    else:
        raise ValueError("y_label must be a string or a list of 3 strings")
    if isinstance(y_value, list):
        y_value_pbmc, y_value_gc, y_value_pb = y_value
    else:
        y_value_pbmc = y_value
        y_value_gc = y_value
        y_value_pb = y_value
    # keep a counter since it won't enumerate when we skip the dummies
    weeks_index = 0
    for ax_index, ax in enumerate(row_axis):
        # remove the dummy axis
        if ax_index == 5 or ax_index == 8:
            ax.remove()
            continue

        if ax_index < 5:
            y_value = y_value_pbmc
        elif ax_index > 5 and ax_index < 8:
            y_value = y_value_gc
        elif ax_index > 8:
            y_value = y_value_pb
        else:
            raise ValueError("ax_index must be between 0 and 9")
        # week lookup is interger of week
        week = timepoints[weeks_index]

        # Plotting df is only the week we are interested in
        plotting_df = plottable_df[plottable_df["weeks_post"] == week]

        # baseline df will take all 0's and replace with lowerbound thresh
        baseline_df = plotting_df.copy()
        baseline_df = baseline_df.replace(0, thresh)
        baseline_df = baseline_df[baseline_df[y_value] == thresh]

        # baseline df may not contain any values for certain timepoints so this just puts in na's if it doesn't have a value for that treatment
        new_dfs = pd.DataFrame(plotting_df["Treatment"].unique(), columns=["Treatment"]).sort_values("Treatment")[::-1]
        baseline_df = pd.concat([baseline_df, new_dfs])

        # replace all 0's with np.nan in the plotting df that we calculate stats on
        if replace_0_with_nan:
            plotting_df = plotting_df.replace(0, math.nan)

        box_df = []
        for _, g_df in plotting_df.groupby(x_value):
            if len(g_df[~g_df[y_value].isna()]) > 3:
                box_df.append(g_df)

        if box_df:
            box_df = pd.concat(box_df).reset_index().sort_values("Treatment")

            # boxplot first
            sns.boxplot(
                ax=ax,
                data=box_df,
                x=x_value,
                y=y_value,
                fliersize=0,
                whis=0,
                palette=palette,
                order=plotting_df["Treatment"].unique(),
            )

        if plot_placebos_seperate:
            placebo_1 = plotting_df[plotting_df["PubID"] == "PubID_001"]
            placebo_2 = plotting_df[plotting_df["PubID"] == "PubID_080"]
            not_placebos = plotting_df[~plotting_df["PubID"].isin(["PubID_001", "PubID_080"])]
        else:
            not_placebos = plotting_df

        # stripplot se

        # make markers circles for not placebos
        sns.stripplot(
            ax=ax,
            data=not_placebos,
            x=x_value,
            y=y_value,
            edgecolor="black",
            linewidth=0.5,
            palette=palette,
            marker="o",
            s=strip_size,
            jitter=value_jitter,
            alpha=value_alpha,
            order=plotting_df["Treatment"].unique(),
        )

        # make markers triangles for placebos
        if plot_placebos_seperate:
            for marker, placebo in [("^", placebo_1), ("s", placebo_2)]:
                sns.stripplot(
                    ax=ax,
                    data=placebo,
                    x=x_value,
                    y=y_value,
                    edgecolor="black",
                    linewidth=0.5,
                    palette=palette,
                    s=strip_size,
                    jitter=value_jitter + jitter_buff,
                    alpha=value_alpha,
                    marker=marker,
                    order=plotting_df["Treatment"].unique(),
                )

        # baseline df is all 0s that are now lower bound thresh
        if plot_placebos_seperate:
            placebo_1 = baseline_df[baseline_df["PubID"] == "PubID_001"]
            placebo_2 = baseline_df[baseline_df["PubID"] == "PubID_080"]
            not_placebos = baseline_df[~baseline_df["PubID"].isin(["PubID_001", "PubID_080"])]
        else:
            not_placebos = baseline_df

        sns.stripplot(
            ax=ax,
            data=not_placebos,
            x=x_value,
            y=y_value,
            edgecolor="black",
            linewidth=0.5,
            jitter=baseline_jitter,
            s=strip_size,
            alpha=baseline_alpha,
            palette=palette,
            order=plotting_df["Treatment"].unique(),
        )
        if plot_placebos_seperate:
            for marker, placebo in [("^", placebo_1), ("s", placebo_2)]:
                sns.stripplot(
                    ax=ax,
                    data=placebo,
                    x=x_value,
                    y=y_value,
                    edgecolor="black",
                    linewidth=0.5,
                    alpha=baseline_alpha,
                    jitter=baseline_jitter,
                    s=strip_size,
                    marker=marker,
                    palette=palette,
                    order=plotting_df["Treatment"].unique(),
                )

        weeks_index += 1

        # writing in manual annotations
        for x in ax.get_xticklabels():
            lookup = x.get_text()
            # get x,y coordinates of the text
            (x1, _) = x.get_position()

            # get the df of just the treatment we are looking at
            treatment_df = plotting_df[plotting_df["Treatment"] == lookup]

            # skip pre-vax
            if week != -4:
                # sucrose will just have a 1:
                if lookup == "DPBS sucrose" and y_text_pos is not None:
                    ax.text(x1, y_text_pos, "1:", fontsize=6, ha="center", fontweight="bold")

                # if everything is not a na
                if not treatment_df[y_value].isna().all() and lookup != "DPBS sucrose" and y_text_pos is not None:
                    median = treatment_df[y_value].median() / 100
                    median = 1 / median
                    median = int(median)
                    if len(str(median)) > 4:
                        budge = -0.1
                    else:
                        budge = 0
                    ax.text(x1 + budge, y_text_pos, median, fontsize=6, ha="center", fontweight="bold")

        ax = adjusts_axis(
            ax, ax_index, week, y_label, annotate_args, plot_x_labels=x_tick_labels, skip_header=skip_header
        )
        ax = adjust_boxplot(ax)
        if scale == "symlog":
            ax.set_yscale("symlog", linthresh=0.00001)
        else:
            ax.set_yscale(scale)
        ax.set_ylim(y_min, y_max)

        def _format(y):
            if y > 1e-4:
                return "{:g}".format(y)
            else:
                if y == 0:
                    return "0"
                if y == 1e-4:
                    return r"$\leq$" + r"$10^{" + str(round(math.log(y, 10))) + r"}$"
                else:
                    return r"$10^{" + str(round(math.log(y, 10))) + r"}$"

        if scale == "symlog" or scale == "log":
            ax.yaxis.set_major_formatter(mtick.FuncFormatter(lambda y, _: _format(y)))
            ax.yaxis.set_major_locator(mtick.LogLocator(base=10, numticks=10))
        if ax_index > 0:
            ax.set_yticklabels([])
        if ax_index == 0 and isinstance(y_label, str):
            ax.set_ylabel(y_label, labelpad=labelpad)
        elif ax_index in [0, 6, 9] and isinstance(y_label, list):
            if ax_index == 0:
                ax.set_ylabel(y_label_pbmc, labelpad=y_label_pad_pbmc)
            elif ax_index == 6:
                ax.set_ylabel(y_label_gc, labelpad=y_label_pad_gc)
            else:
                ax.set_ylabel(y_label_pb, labelpad=y_label_pad_pb)

    return ax


def plot_count_panel(
    row_axis: plt.axes,
    plottable_df: pd.DataFrame,
    combined_timepoints: dict,
    palette: dict,
    annotate_args: dict,
    y_min: int,
    y_max: int,
    x_value: str = "Treatment",
    y_value: str | list[str] = "Number of epitope-specific (KO-GT8++) sequenced IgG BCRs that are VRC01-class",
    y_label: str = "Number of VRC01-class B cells",
    baseline_jitter: int = 0.35,
    baseline_alpha: int = 0.9,
    value_jitter: int = 0.1,
    value_alpha: int = 0.9,
    strip_size: int = 4,
    jitter_buff: int = 0.1,
    scale: str = "symlog",
    skip_header: bool = False,
    plot_x_labels: bool = False,
    major_formattter: mtick.Formatter = None,
    plot_placebos_seperate: bool = True,
    labelpad: int = 4,
):
    """
    Logic that handles sequence counts per timepoint - usually panel A
    """
    # top row panel A is VRC01 class counts
    timepoints = list(combined_timepoints.values())

    if isinstance(y_label, str):
        y_label = y_label
    elif isinstance(y_label, list) and len(y_label) == 3:
        y_label_pbmc, y_label_gc, y_label_pb = y_label
        y_label_pad_pbmc, y_label_pad_gc, y_label_pad_pb = labelpad
    else:
        raise ValueError("y_label must be a string or a list of 3 strings")
    if isinstance(y_value, str):
        y_value = y_value
    elif isinstance(y_value, list) and len(y_label) == 3:
        y_value_pbmc, y_value_gc, y_value_pb = y_value
    else:
        raise ValueError("y_label must be a string or a list of 3 strings")
    # keep a counter since it won't enumerate when we skip the dummies
    weeks_index = 0
    for ax_index, ax in enumerate(row_axis):
        # remove the dummy axis
        if ax_index == 5 or ax_index == 8:
            ax.remove()
            continue

        if (ax_index >= 0 and ax_index <= 5) and isinstance(y_value, list):
            y_value_sub = y_value_pbmc
        elif ax_index in [6, 7] and isinstance(y_value, list):
            y_value_sub = y_value_gc
        elif ax_index == 9 and isinstance(y_value, list):
            y_value_sub = y_value_pb
        else:
            y_value_sub = y_value

        # week lookup is interger of week
        week = timepoints[weeks_index]

        # Plotting df is only the week we are interested in
        plotting_df = plottable_df[plottable_df["weeks_post"] == week]

        box_df = []
        for _, g_df in plotting_df.groupby(x_value):
            if len(g_df[~g_df[y_value_sub].isna()]) > 3:
                box_df.append(g_df)
            else:
                print(f"Skipping, {week},{x_value}")
        if box_df:
            box_df = pd.concat(box_df)
            # boxplot first
            sns.boxplot(
                ax=ax,
                data=box_df,
                x=x_value,
                y=y_value_sub,
                fliersize=0,
                whis=0,
                palette=palette,
                order=plotting_df[x_value].unique(),
            )

        if plot_placebos_seperate:
            placebo_1 = plotting_df[plotting_df["PubID"] == "PubID_001"]
            placebo_2 = plotting_df[plotting_df["PubID"] == "PubID_080"]
            not_placebos = plotting_df[~plotting_df["PubID"].isin(["PubID_001", "PubID_080"])]
        else:
            not_placebos = plotting_df
        # make markers circles for not placebos
        sns.stripplot(
            ax=ax,
            data=not_placebos.query(f"`{y_value_sub}` > 0"),
            x=x_value,
            y=y_value_sub,
            edgecolor="black",
            linewidth=0.5,
            palette=palette,
            marker="o",
            s=strip_size,
            jitter=value_jitter,
            alpha=value_alpha,
            order=plotting_df[x_value].unique(),
        )

        # make markers triangles for placebos
        if plot_placebos_seperate:
            for marker, placebo in zip(["^", "s"], [placebo_1, placebo_2]):
                sns.stripplot(
                    ax=ax,
                    data=placebo.query(f"`{y_value_sub}` > 0"),
                    x=x_value,
                    y=y_value_sub,
                    edgecolor="black",
                    linewidth=0.5,
                    palette=palette,
                    s=strip_size,
                    jitter=value_jitter + jitter_buff,
                    alpha=value_alpha,
                    marker=marker,
                    order=plotting_df[x_value].unique(),
                )

        # plot baseline sep
        sns.stripplot(
            ax=ax,
            data=not_placebos.query(f"`{y_value_sub}` == 0"),
            x=x_value,
            y=y_value_sub,
            edgecolor="black",
            linewidth=0.5,
            palette=palette,
            marker="o",
            s=strip_size,
            jitter=baseline_jitter,
            alpha=baseline_alpha,
            order=plotting_df[x_value].unique(),
        )

        # make markers triangles for placebos baselines
        if plot_placebos_seperate:
            for marker, placebo in zip(["^", "s"], [placebo_1, placebo_2]):
                sns.stripplot(
                    ax=ax,
                    data=placebo.query(f"`{y_value}` == 0"),
                    x=x_value,
                    y=y_value_sub,
                    edgecolor="black",
                    linewidth=0.5,
                    palette=palette,
                    s=strip_size,
                    jitter=baseline_jitter,
                    alpha=baseline_alpha,
                    marker=marker,
                    order=plotting_df["Treatment"].unique(),
                )

        # stripplot se
        weeks_index += 1

        ax = adjusts_axis(
            ax, ax_index, week, y_label, annotate_args, plot_x_labels=plot_x_labels, skip_header=skip_header
        )
        ax = adjust_boxplot(ax)
        # this is a little weird, but essenttially 1-10 is linear > 10 is log, the custom SymLocLocator helps with that
        if scale == "symlog":
            ax.set_ylim(y_min, y_max)
            ax.set_yscale("symlog", linthresh=10)
            if major_formattter:
                ax.yaxis.set_major_formatter(major_formattter)
            else:
                ax.yaxis.set_major_formatter(mtick.ScalarFormatter())  # 10 and 100 are major ticks on a sym log

            ax.yaxis.set_minor_locator(MinorSymLogLocator(10))  # custom class to handle sym log ticks.
        else:
            ax.set_yscale("log")
            ax.set_ylim(y_min, y_max)
        if ax_index > 0:
            ax.set_yticklabels([])
        if ax_index == 0 and isinstance(y_label, str):
            ax.set_ylabel(y_label)
        elif ax_index in [0, 6, 9] and isinstance(y_label, list):
            if ax_index == 0:
                ax.set_ylabel(y_label_pbmc, labelpad=y_label_pad_pbmc)
            elif ax_index == 6:
                ax.set_ylabel(y_label_gc, labelpad=y_label_pad_gc)
            else:
                ax.set_ylabel(y_label_pb, labelpad=y_label_pad_pb)

    return row_axis


def plot_flow_frequencies(data: Data) -> plt.figure:
    # get 3 axes but with dummyies in between in order to add some psuedo spacing
    fig, (row_ax1, row_ax2, row_ax3) = plt.subplots(
        3,
        10,
        figsize=(8.125, 6.8824),
        gridspec_kw={
            "width_ratios": [1, 1, 1, 1, 1, 0.4, 1, 1, 0.4, 1],
        },
        sharex=False,
        sharey=False,
    )

    frequency_df = data.get_flow_and_frequency_data()
    combined_timepoints = data.plot_parameters.get_combined_timepoints()
    gc_timepoints = data.plot_parameters.get_gc_visit_lookup().keys()  # noqa
    frequency_pallete = data.plot_parameters.get_treatment_pallete()
    small_annotate_args = data.plot_parameters.get_small_annotate_args()

    "top row: % of IgG that are GT8++"
    # B cell row - Percent of IgG+ B cells that are GT8++ (without regard to KO binding status)
    # GC cell row - Percent of IgG+ GC B cells that are GT8++ (without regard to KO binding status)
    not_gc_df = frequency_df.query("Visit not in @gc_timepoints")
    gc_df = frequency_df.query("Visit in @gc_timepoints")
    not_gc_df = not_gc_df.rename(
        columns={
            "Percent of IgG+ B cells that are GT8++ (without regard to KO binding status)": "contrived_a",
            "Percent of IgG+ B cells that are epitope-specific (KO-GT8++)": "contrived_b",
            "Percent of GT8++IgG+ B cells that are KO-": "contrived_c",
        }
    ).reset_index(drop=True)
    gc_df = gc_df.rename(
        columns={
            "Percent of IgG+ GC B cells that are GT8++ (without regard to KO binding status)": "contrived_a",
            "Percent of IgG+ GC B cells that are epitope-specific (KO-GT8++)": "contrived_b",
            "Percent of GT8++IgG+ B cells that are KO-": "contrived_c",  # this one is the same but it's gated correctly?
        },
    ).reset_index(drop=True)
    contrived_df = pd.concat([not_gc_df, gc_df])
    plot_frequency_panel(
        row_ax1,
        contrived_df,
        combined_timepoints,
        frequency_pallete,
        small_annotate_args,
        0.00008,
        200,
        "Treatment",
        "contrived_a",
        y_label=[
            "% GT8$^{++}$ among\nIgG$^{+}$ B cells",
            "% GT8$^{++}$ among\nIgG$^{+}$ GC B cells",
            "% GT8$^{++}$ among\nIgD$^{-}$ plasmablasts",
        ],
        y_text_pos=None,
        skip_header=False,
        scale="log",
        thresh=0.0001,
        labelpad=[-4, 0, 0],
        plot_placebos_seperate=False,
    )

    "middle row: % of IgG that are CD4bs-specific"
    plot_frequency_panel(
        row_ax2,
        contrived_df,
        combined_timepoints,
        frequency_pallete,
        small_annotate_args,
        0.00008,
        200,
        "Treatment",
        "contrived_b",
        y_label=[
            "% CD4bs-specific among\nIgG$^{+}$ B cells",
            "% CD4bs-specific among\nIgG$^{+}$ GC B Cells",
            "% CD4bs-specific among\nIgD$^{-}$ plasmablasts",
        ],
        thresh=0.0001,
        y_text_pos=None,
        plot_placebos_seperate=False,
        labelpad=[-4, 0, 0],
    )

    # bottom row: % of GT8++ that are KO-"
    plot_frequency_panel(
        row_ax3,
        contrived_df,
        combined_timepoints,
        frequency_pallete,
        small_annotate_args,
        -0.01,
        109,
        "Treatment",
        "contrived_c",
        y_label=[
            "% KO$^{-}$ among\nGT8$^{++}$IgG$^{+}$ B cells",
            "% KO$^{-}$ among\nGT8$^{++}$IgG$^{+}$ GC B cells",
            "% KO$^{-}$ among\nGT8$^{++}$IgD$^{-}$ plasmablasts",
        ],
        y_text_pos=None,
        x_tick_labels=True,
        thresh=0.01,
        scale="linear",
        labelpad=[10, 0, 0],
        plot_placebos_seperate=False,
    )
    fig.subplots_adjust(hspace=0.1, left=0.11, right=0.975, top=0.9)
    return fig


def plot_count_frequency_response(
    data: Data,
) -> plt.figure:
    """Figure 2 of Paper 1. Plotting FACS counts of epitope specific, frequencies of VRC01 class and vaccine response.
    stratified by timepoint and treatment"""

    # get 3 axes but with dummyies in between in order to add some psuedo spacing
    fig, (row_ax1, row_ax2, row_ax3) = plt.subplots(
        3,
        10,
        figsize=(8.5, 7.2),
        gridspec_kw={
            "width_ratios": [1, 1, 1, 1, 1, 0.4, 1, 1, 0.4, 1],
            "height_ratios": [1, 1, 1],
        },
        sharex=False,
        sharey=False,
    )
    frequency_df = data.get_flow_and_frequency_data()
    plot_params = data.plot_parameters
    combined_timepoints = plot_params.get_combined_timepoints()
    frequency_pallete = plot_params.get_treatment_pallete()
    small_annotate_args = plot_params.get_small_annotate_args()
    gc_timepoints = plot_params.get_gc_visit_lookup().keys()  # noqa
    not_gc_df = frequency_df.query("Visit not in @gc_timepoints")
    gc_df = frequency_df.query("Visit in @gc_timepoints")
    not_gc_df = not_gc_df.rename(
        columns={
            "Percent of IgG+ B cells detected as VRC01-class (missing seq to 0)": "contrived_a",
        }
    ).reset_index(drop=True)
    gc_df = gc_df.rename(
        columns={
            "Percent of IgG+ GC B cells detected as VRC01-class (missing seq to 0)": "contrived_a",
        },
    ).reset_index(drop=True)
    contrived_df = pd.concat([not_gc_df, gc_df])
    # plot pannel A which is sequencing count of epitope specific B cells
    plot_count_panel(
        row_ax1,
        frequency_df,
        combined_timepoints,
        frequency_pallete,
        small_annotate_args,
        -0.5,
        150,
        "Treatment",
        "Number of epitope-specific (KO-GT8++) sequenced IgG BCRs that are VRC01-class",
        y_label=[
            "Number of VRC01-class\nIgG$^{+}$ B cells",
            "Number of VRC01-class\nIgG$^{+}$ GC B cells",
            "Number of VRC01-class\nIgD$^{-}$ plasmablasts",
        ],
        labelpad=[9, 0, 0],
    )
    plot_frequency_panel(
        row_ax2,
        contrived_df,
        combined_timepoints,
        frequency_pallete,
        small_annotate_args,
        0.00008,
        100,
        "Treatment",
        "contrived_a",
        y_label=[
            "% of IgG$^{+}$ memory B cells\ndetected as VRC01-class",
            "% of IgG$^{+}$ GC B cells\ndetected as VRC01-class",
            "% of IgD$^{-}$ plasmablasts\ndetected as VRC01-class",
        ],
        y_text_pos=50,
        labelpad=[-5, 1, 1],
        replace_0_with_nan=True,
    )

    plot_response_panel(
        row_ax3,
        frequency_df,
        combined_timepoints,
        frequency_pallete,
        small_annotate_args,
        0,
        1.05,
        "Treatment",
        "Response (missing seq to 0)",
        y_label=[
            "% of participants with\nVRC01-class BCRs",
            "% of participants with\nVRC01-class BCRs",
            "% of participants with\nVRC01-class BCRs",
        ],
        labelpad=[10, 4, 4],
    )

    fig.subplots_adjust(hspace=0.1, left=0.11, right=0.975, top=0.9)
    return fig


def plot_vrc01_among_epitope_specific(
    data: Data,
) -> plt.figure:
    """Main Figure 4, The percentage of IgG+ B cells that are VRC01-class BCRs

    Parameters
    ----------
    data : Data
    """
    # get 1 axes but with dummyies in between in order to add some psuedo spacing
    fig, (row_ax1, row_ax2) = plt.subplots(
        2,
        10,
        figsize=(8.5, 6.5),
        gridspec_kw={
            "width_ratios": [1, 1, 1, 1, 1, 0.4, 1, 1, 0.4, 1],
        },
        sharex=False,
        sharey=False,
    )

    frequency_df = data.get_flow_and_frequency_data()
    plot_params = data.plot_parameters
    combined_timepoints = plot_params.get_combined_timepoints()
    frequency_pallete = plot_params.get_treatment_pallete()
    small_annotate_args = plot_params.get_small_annotate_args()

    "top row: % of VRC01 class amont eptipope specific"
    plot_frequency_panel(
        row_ax1,
        frequency_df,
        combined_timepoints,
        frequency_pallete,
        small_annotate_args,
        0,
        108,
        "Treatment",
        "Percent of epitope-specific (KO-GT8++) sequenced IgG BCRs that are VRC01-class",
        y_text_pos=None,
        y_label=[
            "% VRC01-class among CD4bs-specific\nIgG$^{+}$ memory BCR sequences",
            "% VRC01-class among CD4bs-specific\nIgG$^{+}$ GC BCR sequences",
            "% VRC01-class among CD4bs-specific\nIgD$^{-}$ PB BCR sequences",
        ],
        skip_header=False,
        thresh=0,
        scale="linear",
        x_tick_labels=False,
        labelpad=[0, 0, 0],
    )
    plot_frequency_panel(
        row_ax2,
        frequency_df,
        combined_timepoints,
        frequency_pallete,
        small_annotate_args,
        0,
        108,
        "Treatment",
        "Percent of GT8++ IgG+ B cells detected as VRC01-class (missing seq to 0)",
        y_text_pos=None,
        y_label=[
            "% VRC01-class among GT8-specific\nIgG$^{+}$ memory BCR sequences",
            "% VRC01-class among GT8-specific\nIgG$^{+}$ GC BCR sequences",
            "% VRC01-class among GT8-specific\nIgD$^{-}$ PB BCR sequences",
        ],
        skip_header=True,
        thresh=0,
        scale="linear",
        x_tick_labels=True,
        labelpad=[0, 0, 0],
    )
    fig.subplots_adjust(hspace=0.1, left=0.11, right=0.97, top=0.85, bottom=0.1)
    return fig


# @todo for v 1.1
# @dataclass
# class FrequencyPlot:
#     row_ax: plt.axes
#     dataframe: pd.DataFrame
#     color_palette: dict[str, str]
#     y_value_pbmc: str
#     y_value_gc: str
#     y_value_pb: str
#     x_value: str
#     hue_value: str
#     skip_header: bool
#     skip_x_labels: bool
#     plot_placebos_seperate: bool
#     y_labels: list[str]
#     y_min: float
#     y_max: float
#     combined_timepoints_and_weeks: dict[str, int]
#     order: list[str]
#     threshold: float = 0.0001

#     def plot(self):
#         # keep a counter since it won't enumerate when we skip the dummies
#         weeks_index = 0

#         # just a list of the weeks in order
#         weeks: list[int] = list(self.combined_timepoints_and_weeks.values())

#         # enumerate through each ax
#         for ax_index, ax in enumerate(self.row_ax):
#             if ax_index == 5 or ax_index == 8:
#                 ax.remove()
#                 continue

#             if ax_index < 5:
#                 y_value = self.y_value_pbmc
#             elif ax_index > 5 and ax_index < 8:
#                 y_value = self.y_value_gc
#             elif ax_index > 8:
#                 y_value = self.y_value_pb
#             else:
#                 raise ValueError("ax_index must be between 0 and 9")

#             # get just the week of interest
#             week = weeks[weeks_index]
#             # Plotting df is only the week we are interested in
#             plotting_df: pd.DataFrame = self.dataframe[self.dataframe["weeks_post"] == week]

#             # baseline df will take all 0's and replace with lowerbound thresh
#             baseline_df: pd.DataFrame = (
#                 plotting_df.replace(0, self.threshold).query(f"`{y_value}` == {self.threshold}").copy()
#             )

#             # baseline df may not contain any values for certain timepoints so this just puts in na's if it doesn't have a value for that treatment
#             all_treatments = pd.DataFrame(plotting_df["Treatment"].unique(), columns=["Treatment"]).sort_values(
#                 "Treatment"
#             )[::-1]
#             baseline_df: pd.DataFrame = pd.concat([all_treatments, baseline_df])

#             # box df will just contain those groups which we will plot
#             box_df: list[pd.DataFrame] = []

#             # only plot things in boxplot that have more than three values
#             for _, g_df in plotting_df.groupby(self.x_value):

#                 if len(g_df[~g_df[y_value].isna()]) > 3:
#                     box_df.append(g_df)

#             # if we have stuff in boxplot, atttempt to plot it
#             if box_df:
#                 box_df = pd.concat(box_df).reset_index().sort_values("Treatment")

#                 # boxplot first
#                 sns.boxplot(
#                     ax=ax,
#                     data=box_df,
#                     x=self.x_value,
#                     y=y_value,
#                     fliersize=0,
#                     whis=0,
#                     palette=self.color_palette,
#                     order=self.order,
#                 )


# def plot_flow_frequencies(data: Data) -> plt.figure:
#     """Figure 1. Plotting FACS frequencies stratified by timepoint and treatment"""

#     # get 3 axes but with dummyies in between in order to add some psuedo spacing
#     _width_ratios = [1, 1, 1, 1, 1, 0.4, 1, 1, 0.4, 1]
#     _fig_size = (8.125, 6.8824)

#     fig, (row_ax1, row_ax2, row_ax3) = plt.subplots(
#         3,
#         10,
#         figsize=_fig_size,
#         gridspec_kw={"width_ratios": _width_ratios},
#         sharex=False,
#         sharey=False,
#     )

#     # frequency and sequence dataframe
#     frequency_df = data.get_flow_and_frequency_data()

#     # get vist and weeks of everything
#     combined_timepoints_and_weeks = data.plot_parameters.get_combined_timepoints()

#     # just get the visit ids of the gc
#     gc_timepoints = data.plot_parameters.get_gc_visit_lookup().keys()  # noqa

#     # get the palette we use
#     frequency_pallete = data.plot_parameters.get_treatment_pallete()

#     # some frequent annotation arguments
#     small_annotate_args = data.plot_parameters.get_small_annotate_args()

#     # First, break up gc and non-gc
#     not_gc_df = frequency_df.query("Visit not in @gc_timepoints")
#     gc_df = frequency_df.query("Visit in @gc_timepoints")

#     # rename the columns so we can plot them using the same y-value (contrived_a,b,c)
#     not_gc_df = not_gc_df.rename(
#         columns={
#             "Percent of IgG+ B cells that are GT8++ (without regard to KO binding status)": "contrived_a",
#             "Percent of IgG+ B cells that are epitope-specific (KO-GT8++)": "contrived_b",
#             "Percent of GT8++IgG+ B cells that are KO-": "contrived_c",
#         }
#     ).reset_index(drop=True)
#     gc_df = gc_df.rename(
#         columns={
#             "Percent of IgG+ GC B cells that are GT8++ (without regard to KO binding status)": "contrived_a",
#             "Percent of IgG+ GC B cells that are epitope-specific (KO-GT8++)": "contrived_b",
#             "Percent of GT8++IgG+ B cells that are KO-": "contrived_c",  # this one is the same but it's gated correctly?
#         },
#     ).reset_index(drop=True)

#     # put them back together
#     contrived_df = pd.concat([not_gc_df, gc_df])

#     # Top Plot, IgG Antigen Spcific (GT8++)
#     frequency_plot_top_row = FrequencyPlot(
#         row_ax=row_ax1,
#         dataframe=frequency_df,
#         color_palette=frequency_pallete,
#         y_value_pbmc="Percent of IgG+ B cells that are GT8++ (without regard to KO binding status)",
#         y_value_pb="Percent of IgG+ B cells that are GT8++ (without regard to KO binding status)",
#         y_value_gc="Percent of IgG+ GC B cells that are epitope-specific (KO-GT8++)",
#         x_value="Treatment",
#         hue_value="Treatment",
#         y_labels=[
#             "% GT8$^{++}$ among\nIgG$^{+}$ B cells",
#             "% GT8$^{++}$ among\nIgG$^{+}$ GC B cells",
#             "% GT8$^{++}$ among\nIgD$^{-}$ plasmablasts",
#         ],
#         skip_header=False,
#         skip_x_labels=False,
#         plot_placebos_seperate=False,
#         y_min=0.00008,
#         y_max=200,
#         combined_timepoints_and_weeks=combined_timepoints_and_weeks,
#         order=["DPBS sucrose", "20 µg eOD-GT8 60mer + AS01B", "100 µg eOD-GT8 60mer + AS01B"],
#     )
#     frequency_plot_top_row.plot()
#     return fig

#     # "top row:O % of IgG that are GT8++ (Jimmy PT figs 11, 12)"
#     # # B cell row - Percent of IgG+ B cells that are GT8++ (without regard to KO binding status)
#     # # GC cell row - Percent of IgG+ GC B cells that are GT8++ (without regard to KO binding status)
#     # contrived_df = pd.concat([not_gc_df, gc_df])
#     # _plot_frequency_panel(
#     #     row_ax1,
#     #     contrived_df,
#     #     combined_timepoints,
#     #     frequency_pallete,
#     #     small_annotate_args,
#     #     0.00008,
#     #     200,
#     #     "Treatment",
#     #     "contrived_a",
#     #     y_label=[
#     #         "% GT8$^{++}$ among\nIgG$^{+}$ B cells",
#     #         "% GT8$^{++}$ among\nIgG$^{+}$ GC B cells",
#     #         "% GT8$^{++}$ among\nIgD$^{-}$ plasmablasts",
#     #     ],
#     #     y_text_pos=None,
#     #     skip_header=False,
#     #     scale="log",
#     #     thresh=0.0001,
#     #     labelpad=[-4, 0, 0],
#     #     plot_placebos_seperate=False,
#     # )

#     # "middle row: % of IgG that are CD4bs-specific  (Jimmy PT figs 7, 8)"
#     # _plot_frequency_panel(
#     #     row_ax2,
#     #     contrived_df,
#     #     combined_timepoints,
#     #     frequency_pallete,
#     #     small_annotate_args,
#     #     0.00008,
#     #     200,
#     #     "Treatment",
#     #     "contrived_b",
#     #     y_label=[
#     #         "% CD4bs-specific among\nIgG$^{+}$ B cells",
#     #         "% CD4bs-specific among\nIgG$^{+}$ GC B Cells",
#     #         "% CD4bs-specific among\nIgD$^{-}$ plasmablasts",
#     #     ],
#     #     thresh=0.0001,
#     #     y_text_pos=None,
#     #     plot_placebos_seperate=False,
#     #     labelpad=[-4, 0, 0],
#     # )

#     # "bottom row: % of GT8++ that are KO- (Jimmy PT figs 13, 14)"
#     # "columna: Percent of GT8++IgG+ B cells that are KO-"
#     # _plot_frequency_panel(
#     #     row_ax3,
#     #     contrived_df,
#     #     combined_timepoints,
#     #     frequency_pallete,
#     #     small_annotate_args,
#     #     -0.01,
#     #     109,
#     #     "Treatment",
#     #     "contrived_c",
#     #     y_label=[
#     #         "% KO$^{-}$ among\nGT8$^{++}$IgG$^{+}$ B cells",
#     #         "% KO$^{-}$ among\nGT8$^{++}$IgG$^{+}$ GC B cells",
#     #         "% KO$^{-}$ among\nGT8$^{++}$IgD$^{-}$ plasmablasts",
#     #     ],
#     #     y_text_pos=None,
#     #     x_tick_labels=True,
#     #     thresh=0.01,
#     #     scale="linear",
#     #     labelpad=[10, 0, 0],
#     #     plot_placebos_seperate=False,
#     # )
#     # plt.subplots_adjust(hspace=0.1, left=0.11, right=0.975, top=0.9)
#     # plt.savefig(outpath + ".pdf")
#     # plt.savefig(outpath + ".png", dpi=300)
#     # plt.savefig(outpath + ".svg")
