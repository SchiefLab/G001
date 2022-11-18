import math
import re
from itertools import product
import pandas as pd
import seaborn as sns
from pathlib import Path
from Levenshtein._levenshtein import distance
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import ticker as mtick
from matplotlib.patches import Patch
from sadie.airr.airrtable import LinkedAirrTable

from g001.figures.somatic import format_xtick_labels
from g001.figures.util import adjust_boxplot
from g001.data import Data

from g001.figures.routines import add_mutational_sets
from g001.figures.seq_logo import plot_seq_logo_panel, plot_sequence_logo

explicit_mapping = {
    "IGKV1-33": "#91FCC0",
    "IGKV1-5": "#E377C2",
    "IGKV3-20": "#2078B4",
    "IGLV2-14": "#17BFD0",
    "IGKV3-15": "#9567BD",
    "IGLV2-11": "#BDBD23",
    "IGLV2-23": "#FF7F0F",
    "Other": "#F2F0F2",
}

sort_order = [
    "IGKV1-33",
    "IGKV3-20",
    "IGKV1-5",
    "IGKV3-15",
    "IGLV2-14",
    "IGLV2-23",
    "IGLV2-11",
    "Other",
]


def plot_sequence_features(data: Data) -> plt.figure:
    sequence_df = data.get_unblinded_sequences()
    plot_parm = data.plot_parameters
    dekosky = data.get_dekosky()
    kappa_vrc01_class = data.get_kappa_vrc01_aa()
    lambda_vrc01_class = data.get_lambda_vrc01_aa()
    oas_5_len = data.get_oas_five_len()
    oas_vh12 = data.get_oas_vh12()
    vrc01_ref_airr = data.get_vr01_ref_airr()
    human_naive_5_len = data.get_human_naive_airr().query("is_vrc01_class")
    palette = plot_parm.get_pallete()

    min_set = [
        "12A12",
        "12A21",
        "N6",
        "VRC27",
        "N17",
        "N60P1.1",
        "N60P25.1",
        "N60P23",
        "PCIN63_71I",
        "PCIN63_66B",
        "PCIN63_71G",
        "NIH45-46",
        "VRC07b",
        "VRC23",
        "VRC01",
        "VRC02",
        "VRC18",
        "VRC08",
        "VRC-PG19",
    ]

    # establish A4 size
    fig = plt.figure(figsize=(8.5, 11.0))

    plt.rcParams["axes.linewidth"] = 0.5
    plt.rcParams["xtick.minor.width"] = 0.5
    plt.rcParams["ytick.minor.width"] = 0.5
    plt.rcParams["xtick.major.width"] = 0.5
    plt.rcParams["ytick.major.width"] = 0.5

    # panel start left margin
    panel_top = 0.95
    panels_left = 0.1
    panels_right = 0.98

    # how high should an indivual panel b
    indidual_panel_height = 0.12

    # where will seq logo start?
    seq_logo_panel_left = 0.70
    left_right_buffer = 0.08

    # multi panel start is 5 panels high
    multi_panel_bottom = panel_top - (5 * indidual_panel_height)

    # panel A,D,E,F
    multi_panel_grid_spec = fig.add_gridspec(
        ncols=4,  # combiend vrc01, combined nonvrc01, control, legend
        nrows=4,
        top=panel_top,
        left=panels_left,
        right=seq_logo_panel_left - left_right_buffer,
        bottom=multi_panel_bottom,
        wspace=0.09,
        hspace=0.75,
        width_ratios=[1, 1, 1.7 / 5, 1.2 / 5],
    )

    plottable_df, plottable_with_muts_df, combined_5_len_oas_df, oas_vh12_df = add_sequence_metrics(
        sequence_df, oas_5_len, vrc01_ref_airr, oas_vh12
    )

    # plot light chain usage
    plot_light_segment_usage(fig, multi_panel_grid_spec, plottable_df, dekosky, explicit_mapping)

    # use Q/E ?
    plot_usage_metric(
        fig, multi_panel_grid_spec, plottable_df, oas_vh12_df, combined_5_len_oas_df, palette, metric="has_QE_at_96"
    )

    # dist to lc
    plot_light_chain_dist(fig, multi_panel_grid_spec, plottable_df, combined_5_len_oas_df, palette)

    # use W100b?
    plot_usage_metric(
        fig, multi_panel_grid_spec, plottable_df, oas_vh12_df, combined_5_len_oas_df, palette, metric="has_w100b"
    )

    # Seq logo panel B,C
    seq_logo_grid_spec = fig.add_gridspec(
        11,
        2,
        left=seq_logo_panel_left,
        right=panels_right,
        bottom=multi_panel_bottom,
        top=panel_top + 0.015,
        wspace=0.08,
        height_ratios=[0.3, 1, 1, 1, 1, 0.8, 0.3, 1, 1, 1, 1],
    )
    plot_seq_logo_panel(
        fig,
        seq_logo_grid_spec,
        kappa_vrc01_class,
        lambda_vrc01_class,
        combined_5_len_oas_df,
        plottable_df,
        human_naive_5_len,
    )

    # panel e  is cloud plot
    panel_widths = [1, 1, 1.5 / 6]
    panel_e_gs = fig.add_gridspec(
        1,
        3,
        width_ratios=panel_widths,
        left=panels_left,
        right=panels_right,
        top=multi_panel_bottom - 0.07,
        bottom=multi_panel_bottom - 0.07 - indidual_panel_height,
        wspace=0.13,
    )

    # cloud plot
    plot_key_mutations(
        fig,
        panel_e_gs,
        palette,
        plottable_with_muts_df,
        vrc01_ref_airr[vrc01_ref_airr["sequence_id"].isin(min_set)].copy(),
        data.get_cottrell_path(),
    )
    return fig


def manipulate_df(
    df,
    key=[
        "pubid",
        "timepoint",
        "dose_group",
        "weeks_post",
        "locus_light",
        "is_vrc01_class",
        "vaccine_group",
    ],
):
    df["distance_to_known_lcdr3"] = df[["distance_to_known_lcdr3_kappa", "distance_to_known_lcdr3_lambda"]].min(axis=1)
    _df = (
        df.groupby(key)
        .apply(report_lcdr3)
        .stack()
        .reset_index()
        .rename(
            {0: "frequency"},
            axis=1,
        )
    )
    return _df


def report_lcdr3(df):
    b = df["distance_to_known_lcdr3"].value_counts().sort_index()
    c = (b / b.sum()).cumsum()
    c = c.reindex([float(i) for i in range(6)]).ffill().fillna(0.0).transpose()
    return c


def find_lowest_len(x, ref):
    return min(list(map(lambda y: distance(x, y), ref)))


def add_sequence_metrics(
    unblinded_df: pd.DataFrame,
    oas_5_len_df: pd.DataFrame,
    vrc01_ref_airr_df: pd.DataFrame,
    oas_vh12_df: pd.DataFrame,
):

    # unblindeed dataframe of all sequences
    plottable_df = unblinded_df.copy()
    plottable_df_with_muts = unblinded_df.copy()

    # add metrics to the dataframe
    # only plot the first call and the gene, not the allele
    plottable_df["light_plot"] = plottable_df["v_call_light"].str.split(",").str.get(0).str.split("*").str.get(0)

    # change 3D-15 to 3-15
    plottable_df.loc[plottable_df[plottable_df["light_plot"] == "IGKV3D-15"].index, "light_plot"] = "IGKV3-15"

    # anything not in our list is "other"
    plottable_df.loc[
        plottable_df[~plottable_df["light_plot"].isin(explicit_mapping.keys())].index,
        "light_plot",
    ] = "Other"

    # add which shot they are at, prime or boost
    plottable_df.loc[plottable_df[plottable_df["weeks_post"] <= 8].index, "shot"] = "first shot"
    plottable_df.loc[plottable_df[plottable_df["weeks_post"] > 8].index, "shot"] = "second shot"

    plottable_df["has_QE_at_96"] = plottable_df["cdr3_aa_light"].apply(
        lambda x: True if len(x) == 5 and x[-2] in ["Q", "E"] else False
    )

    # our data
    plottable_df["has_w100b"] = plottable_df["junction_aa_heavy"].str[-6] == "W"

    # briney sotosoto_briney_vh12["has_w100b"] = soto_briney_vh12["cdr3_aa"].str[-5] == "W"
    oas_vh12_df["has_w100b"] = oas_vh12_df["cdr3_aa"].str[-5] == "W"

    # find distances from ref cdrs
    # distance to nearest lcdr3
    vrc01_ref_airr_kappa_seqs = vrc01_ref_airr_df.query("locus_light == 'IGK'")["cdr3_aa_light"].to_list()
    vrc01_ref_airr_lambda_seqs = vrc01_ref_airr_df.query("locus_light == 'IGL'")["cdr3_aa_light"].to_list()

    # Find all distances to nearest lcdr3 for oas
    oas_5_left_ref_airr = oas_5_len_df.query("oas_subject != 'None'")
    oas_5_left_ref_airr_kappa = oas_5_left_ref_airr.query("locus == 'IGK'").copy()
    oas_5_left_ref_airr_lambda = oas_5_left_ref_airr.query("locus == 'IGL'").copy()
    oas_5_left_ref_airr_kappa["distance_to_known_lcdr3_kappa"] = oas_5_left_ref_airr_kappa["cdr3_aa"].apply(
        lambda x: find_lowest_len(x, vrc01_ref_airr_kappa_seqs)
    )
    oas_5_left_ref_airr_lambda["distance_to_known_lcdr3_lambda"] = oas_5_left_ref_airr_lambda["cdr3_aa"].apply(
        lambda x: find_lowest_len(x, vrc01_ref_airr_lambda_seqs)
    )

    # combine them backtogeter
    combined_oas_5_len_df = pd.concat([oas_5_left_ref_airr_kappa, oas_5_left_ref_airr_lambda]).reset_index(drop=True)
    combined_oas_5_len_df["has_QE_at_96"] = combined_oas_5_len_df["cdr3_aa"].apply(
        lambda x: True if x[-2] in ["Q", "E"] else False
    )

    # so far only probe specific
    class_df = plottable_df.query("weeks_post>0").query("vaccine_group == 'vaccine'")
    class_df_kappa = class_df.query("locus_light == 'IGK'").copy()
    class_df_lambda = class_df.query("locus_light == 'IGL'").copy()
    class_df_kappa["distance_to_known_lcdr3_kappa"] = class_df_kappa["cdr3_aa_light"].apply(
        lambda x: find_lowest_len(x, vrc01_ref_airr_kappa_seqs)
    )
    class_df_lambda["distance_to_known_lcdr3_lambda"] = class_df_lambda["cdr3_aa_light"].apply(
        lambda x: find_lowest_len(x, vrc01_ref_airr_lambda_seqs)
    )
    combined_light_df = pd.concat([class_df_kappa, class_df_lambda]).reset_index(drop=True)

    # return unblinded df, mutational_df, all oas df with 5 len and combined light_df with lookup
    return combined_light_df, plottable_df_with_muts, combined_oas_5_len_df, oas_vh12_df


def get_top_n_percet(df, percent=0.9):
    """Groupby funciton to get qunatile of key residues"""
    return pd.Series({"residues": df["cottrell_focused_v_common_score"].quantile(percent, interpolation="midpoint")})


def process_dekosky(dekosky_df: pd.DataFrame, explicit_mapping: dict) -> pd.DataFrame:
    # add dekosky data
    dekosky_df.loc[
        dekosky_df[dekosky_df["v_call_top_light"] == "IGKV3D-15"].index,
        "v_call_top_light",
    ] = "IGKV3-15"
    dekosky_df.loc[
        dekosky_df[~dekosky_df["v_call_top_light"].isin(explicit_mapping.keys())].index,
        "v_call_top_light",
    ] = "Other"

    deskosky_ready = (
        dekosky_df.groupby(["Replicate", "donor"])
        .apply(lambda x: x["v_call_top_light"].value_counts(normalize=True))
        .reset_index()
        .rename({"v_call_top_light": "light_plot"}, axis=1)
        .groupby(["level_2"])
        .mean()
        .drop("Replicate", axis=1)
        .reset_index()
    )
    deskosky_ready["x-axis"] = "Dekosky"
    return deskosky_ready


def plot_light_segment_usage(
    fig: plt.figure, panel_a_gs: plt.GridSpec, plottable_df: pd.DataFrame, dekosky: pd.DataFrame, explicit_mapping: dict
):
    """Panel A is light chain usage among VRC01c and nonVRC01c along with the dekosky control"""

    # make three subplots
    light_usage_post_vaccine_vrc01 = fig.add_subplot(panel_a_gs[0, 0])
    light_usage_post_vaccine_nonvrc01 = fig.add_subplot(panel_a_gs[0, 1])
    dekosky_axis = fig.add_subplot(panel_a_gs[0, 2])
    legend_axis = fig.add_subplot(panel_a_gs[0, 3], frameon=False)

    def _pivot(df: pd.DataFrame, weeks_in: list[int], explicit_mapping: dict[str, str]):
        # multiindex is all weeks * all possible light chains,
        mi = pd.MultiIndex.from_product([weeks_in, explicit_mapping.keys()], names=["weeks_post", "light_chain"])

        # pivot and reindex. the mi will take care of the values we don't find and fill in an NA
        pivot_df = (
            df.groupby(["weeks_post"])
            .apply(lambda x: x["light_plot"].value_counts(normalize=True))
            .to_frame()
            .reindex(mi)
            .rename({"light_plot": "freq"}, axis=1)
            .reset_index()
            .pivot("weeks_post", "light_chain", "freq")
        )
        print(pivot_df)
        return pivot_df

    # always reset index, you never know
    plottable = plottable_df.copy().reset_index(drop=True)

    # throw out pre vaccine and placebo
    post_vaccine = plottable[(plottable["timepoint"] != "V02") & (plottable["vaccine_group"] == "vaccine")]
    post_vaccine.loc[post_vaccine[post_vaccine["weeks_post"] <= 8].index, "shot"] = "first shot"
    post_vaccine.loc[post_vaccine[post_vaccine["weeks_post"] > 8].index, "shot"] = "second shot"

    low_dose_label = r"20 $\mu$g"
    high_dose_label = r"100 $\mu$g"

    # combined vrc01 + nonvrc01
    combined = (
        post_vaccine.groupby(["dose_group", "is_vrc01_class", "shot", "has_vh1-2"])
        .apply(lambda x: x["light_plot"].value_counts(normalize=True))
        .to_frame()
        .reset_index()
        .sort_values(["is_vrc01_class", "dose_group", "shot"], ascending=[0, 0, 0])
    )

    # split into vrc01 and nonvrc01
    combined_vrc01 = combined.query("is_vrc01_class == True")
    combined_nonvrc01 = combined.query("is_vrc01_class != True").query("`has_vh1-2`")

    # data dump
    # combined_vrc01.reset_index(drop=True).to_feather("./notebooks/data/combined_vrc01_light_usage.feather")
    # combined_nonvrc01.reset_index(drop=True).to_feather("./notebooks/data/combined_nonvrc01_light_usage.feather")
    for ax, df in [
        (light_usage_post_vaccine_vrc01, combined_vrc01),
        (light_usage_post_vaccine_nonvrc01, combined_nonvrc01),
    ]:
        pivoted_plot = (
            df.rename({"level_4": "gene"}, axis=1)
            .pivot(["dose_group", "shot"], "gene", "light_plot")
            .fillna(0.0)
            .sort_index(level=[0], ascending=[0])
        )
        colors = [explicit_mapping[i] for i in sort_order]
        pivoted_plot.loc[:, sort_order].plot(
            kind="bar",
            stacked=True,
            color=colors,
            width=0.75,
            linewidth=1,
            edgecolor="black",
            ax=ax,
        )

    # pivot deksoky
    colors = [explicit_mapping[i] for i in sort_order]
    dekosky_processed = process_dekosky(dekosky, explicit_mapping)
    dekosky_pivot = dekosky_processed.pivot("x-axis", "level_2", "light_plot").loc[:, sort_order]
    dekosky_pivot.plot(
        kind="bar",
        stacked=True,
        color=colors,
        width=0.75,
        linewidth=1,
        edgecolor="black",
        ax=dekosky_axis,
    )
    for index, (ax, title) in enumerate(
        [
            (light_usage_post_vaccine_vrc01, "VRC01 Class"),
            (light_usage_post_vaccine_nonvrc01, "Non-VRC01 Class\nwith VH1-2"),
            (dekosky_axis, "Control"),
        ]
    ):
        if index in [0, 1]:
            ax.set_xticklabels(
                [
                    low_dose_label + "\nPrime",
                    low_dose_label + "\nBoost",
                    high_dose_label + "\nPrime",
                    high_dose_label + "\nBoost",
                ],
                rotation=0,
                fontsize=7,
            )
        else:
            ax.set_xticklabels(["DeKosky\nVH1-2"], rotation=0, fontsize=7)
        ax.legend_.remove()
        ax.set_xlabel("")
        if index != 0:
            ax.set_yticklabels([])
            ax.set_ylabel("")
        else:
            label = "% BCRs with VRC01-class \n" + r"bnAb $\mathregular{V_{K/L}}$"
            ax.set_ylabel(label, fontsize=8)
            ax.tick_params(axis="y", labelsize=8)

        ax.yaxis.set_major_locator(mtick.MultipleLocator(0.2))
        ax.yaxis.set_minor_locator(mtick.MultipleLocator(0.1))
        if index == 0:
            ax.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1, symbol=None))
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_ylim(0, 1.02)
        if len(title.split("\n")) > 1:
            ax.set_title(title)
        else:
            ax.set_title(title, pad=10)

    legend_axis.set_xticks([])
    legend_axis.set_yticks([])
    custom_lines = []
    for label in sort_order:
        color = explicit_mapping[label]
        if label != "Other":
            crude_label = label[2:].replace("V", "")
        else:
            crude_label = label
        custom_lines.append(Patch(facecolor=color, edgecolor="black", linewidth=1, label=crude_label)),
    legend_axis.legend(
        custom_lines,
        [i._label for i in custom_lines],
        loc="upper center",
        frameon=False,
        handlelength=0.8,
        ncol=1,
        fontsize=8,
        bbox_to_anchor=(0.5, 1.05),
        labelspacing=0.1,
    )


def plot_usage_metric(
    fig: plt.figure,
    panel_a_gs: plt.GridSpec,
    plottable_df: pd.DataFrame,
    oas_vh12_df: pd.DataFrame,
    oas_5_len_df: pd.DataFrame,
    palette: dict,
    metric: str,
):
    if metric == "has_QE_at_96":
        row = 1
    elif metric == "has_w100b":
        row = 3
    else:
        raise ValueError("Invalid metric")

    # make three subplots
    post_vaccine_vrc01 = fig.add_subplot(panel_a_gs[row, 0])
    post_vaccine_nonvrc01 = fig.add_subplot(panel_a_gs[row, 1])
    oas_axis = fig.add_subplot(panel_a_gs[row, 2])
    legend_axis = fig.add_subplot(panel_a_gs[row, 3], frameon=False)

    # always reset index, you never know
    plottable = plottable_df.copy().reset_index(drop=True)

    # throw out pre vaccine and placebo
    post_vaccine = (
        plottable[(plottable["timepoint"] != "V02") & (plottable["vaccine_group"] == "vaccine")]
        .copy()
        .reset_index(drop=True)
    )

    low_dose_label = r"20 $\mu$g"
    high_dose_label = r"100 $\mu$g"

    # combined vrc01 + nonvrc01
    if metric == "has_QE_at_96":
        non_vrc01_metric = "has_5_len_lcdr3"
    else:
        non_vrc01_metric = "has_vh1-2"
    post_vaccine = post_vaccine[post_vaccine[non_vrc01_metric]]

    def _determine_metric(df: pd.DataFrame) -> pd.Series:
        return pd.Series({True: df[metric].sum() / len(df), False: 1 - (df[metric].sum() / len(df))})

    combined = (
        post_vaccine.groupby(["dose_group", "pubid", "is_vrc01_class", "shot"])
        .apply(_determine_metric)
        .stack()
        .reset_index()
        .rename({"level_4": metric, 0: "freq"}, axis=1)
    )

    # if metric == "has_QE_at_96":
    #     combined.reset_index(drop=True).to_feather("./notebooks/data/combined_with_metrics_QE96.feather")
    # else:
    #     combined.reset_index(drop=True).to_feather("./notebooks/data/combined_with_metrics_W100B.feather")

    # split into vrc01 and nonvrc01
    # true is frequency of QE
    # Dose_Group  is_vrc01_class         shot  level_3  has_QE_at_96
    # Low Dose            True  second shot     True      0.613240
    # Low Dose            True  second shot    False      0.386760
    combined["x_axis"] = combined["dose_group"] + "_" + combined["shot"]
    combined_vrc01 = combined.query("is_vrc01_class == True").query(metric)
    combined_nonvrc01 = combined.query("is_vrc01_class != True").query(metric)
    for vrc01_class in [True, False]:
        x_nest = ["Low Dose_first shot", "Low Dose_second shot", "High Dose_first shot", "High Dose_second shot"]
        for shot in ["first shot", "second shot"]:
            for dose in ["Low Dose", "High Dose"]:
                x_axis = dose + "_" + shot
                if vrc01_class:
                    ax = post_vaccine_vrc01
                    c = "vrc01_class"
                    df = combined_vrc01
                else:
                    ax = post_vaccine_nonvrc01
                    c = "nonvrc01_class"
                    df = combined_nonvrc01
                color = palette[dose][c]
                sns.boxplot(
                    data=df.query(f"x_axis == '{x_axis}'"),
                    x="x_axis",
                    y="freq",
                    linewidth=1,
                    color=color,
                    ax=ax,
                    whis=[10, 90],
                    order=x_nest,
                    fliersize=0,
                )
                sns.stripplot(
                    data=df.query(f"x_axis == '{x_axis}'"),
                    x="x_axis",
                    y="freq",
                    linewidth=1,
                    color=color,
                    ax=ax,
                    order=x_nest,
                    size=4,
                    edgecolor="black",
                )
        adjust_boxplot(ax)

    if metric == "has_QE_at_96":

        control_df = (
            oas_5_len_df.groupby(["oas_subject", "oas_author"])
            .apply(lambda x: x["has_QE_at_96"].value_counts(normalize=True))
            .to_frame()
            .reset_index()
            .assign(dose_group="OAS 5 Len")
            .query("level_2")
        )
    else:
        control_df = (
            oas_vh12_df.groupby(["oas_subject", "oas_author"])
            .apply(lambda x: x["has_w100b"].value_counts(normalize=True))[True]
            .to_frame()
            .assign(dose_group="OAS 5 Len")
            .rename({True: metric}, axis=1)
        )
    sns.boxplot(
        data=control_df,
        x="dose_group",
        y=metric,
        linewidth=1,
        color="#F2F0F2",
        ax=oas_axis,
        whis=[10, 90],
        fliersize=0,
    )
    sns.stripplot(
        data=control_df,
        x="dose_group",
        y=metric,
        linewidth=1,
        color="#F2F0F2",
        ax=oas_axis,
        edgecolor="black",
        size=4,
        jitter=0.1,
    )
    adjust_boxplot(oas_axis)
    if metric == "has_QE_at_96":
        y_label = "% BCRs with E/Q at\nLCDR3 position 96"
        non_vrc01_label = "Non-VRC01 Class\nwith 5aa LCDR3"
        control_label = "OAS 5aa L3"
    else:
        y_label = "% BCRs with\n" + r"$\mathregular{Trp_{105-3}}$ in HCDR3"
        # y_label = "% BCRs with W at position -5\nfrom end of HCDR3"
        non_vrc01_label = "Non-VRC01 Class\nwith VH1-2"
        control_label = "OAS VH1-2"
    for index, (ax, title) in enumerate(
        [
            (post_vaccine_vrc01, "VRC01 Class"),
            (post_vaccine_nonvrc01, non_vrc01_label),
            (oas_axis, "Control"),
        ]
    ):
        x_ticks = [
            low_dose_label + "\nPrime",
            low_dose_label + "\nBoost",
            high_dose_label + "\nPrime",
            high_dose_label + "\nBoost",
        ]
        if index in [0, 1]:
            ax.set_xticklabels(
                x_ticks,
                rotation=0,
                fontsize=7,
            )
        else:
            ax.set_xticklabels([control_label], rotation=0, fontsize=7)
        if ax.legend_:
            ax.legend_.remove()
        ax.set_xlabel("")
        if index != 0:
            ax.set_yticklabels([])
            ax.set_ylabel("")
        else:
            ax.set_ylabel(y_label, fontsize=8)
            ax.tick_params(axis="y", labelsize=8)

        ax.yaxis.set_major_locator(mtick.MultipleLocator(0.2))
        ax.yaxis.set_minor_locator(mtick.MultipleLocator(0.1))
        if index == 0:
            ax.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1, symbol=None))
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_ylim(0, 1.05)
        ax.set_title(title)

    legend_axis.axis("off")


def plot_light_chain_dist(
    fig,
    gs,
    combined_df: pd.DataFrame,
    combined_oas_df: pd.DataFrame,
    palette: dict,
):
    """Distances from known LCDRS"""
    thresh = 0

    panel_c_gs = gs
    distance_vrc01_combined = fig.add_subplot(panel_c_gs[2, 0])
    distance_nonvrc01_combined = fig.add_subplot(panel_c_gs[2, 1])
    distance_control = fig.add_subplot(panel_c_gs[2, 2])
    distance_legend = fig.add_subplot(panel_c_gs[2, 3])

    key = [
        "pubid",
        "dose_group",
        "locus_light",
        "is_vrc01_class",
        "vaccine_group",
    ]
    class_all_timepoint = manipulate_df(combined_df.query("has_5_len_lcdr3").copy(), key=key)
    class_all_timepoint["x-axis"] = class_all_timepoint["dose_group"] + "_" + class_all_timepoint["locus_light"]

    # class_all_timepoint.reset_index(drop=True).to_feather("./notebooks/data/dist_to_lcdr3.feather")
    x_axis_order = ["Low Dose_IGK", "Low Dose_IGL", "High Dose_IGK", "High Dose_IGL"]
    for is_igl in [True, False]:
        if is_igl:
            ax = distance_vrc01_combined
            color_lookup = "vrc01_class_kappa"
        else:
            ax = distance_nonvrc01_combined
            color_lookup = "nonvrc01_class_kappa"

        for high_low in ["Low Dose", "High Dose"]:
            for kappa_lambda in ["IGK", "IGL"]:
                plot_df = class_all_timepoint.query(
                    f"locus_light == '{kappa_lambda}' and dose_group == '{high_low}' and is_vrc01_class =={is_igl} and distance_to_known_lcdr3 == {thresh}"
                )
                if high_low == "Low Dose":
                    pal = palette["Low Dose"][color_lookup]
                else:
                    pal = palette["High Dose"][color_lookup]
                if kappa_lambda == "IGK":
                    marker = "o"
                else:
                    marker = "^"
                if not plot_df.query("frequency > 0").empty:
                    sns.boxplot(
                        data=plot_df.query("frequency > 0"),
                        x="x-axis",
                        y="frequency",
                        dodge=True,
                        linewidth=1,
                        color=pal,
                        ax=ax,
                        whis=[10, 90],
                        order=x_axis_order,
                        fliersize=0,
                    )
                    sns.stripplot(
                        data=plot_df.query("frequency > 0"),
                        x="x-axis",
                        y="frequency",
                        dodge=True,
                        edgecolor="black",
                        linewidth=1,
                        color=pal,
                        order=x_axis_order,
                        ax=ax,
                        size=4,
                        marker=marker,
                    )
                    adjust_boxplot(ax)
                if not plot_df.query("frequency  == 0").empty:
                    sns.stripplot(
                        data=plot_df.query("frequency == 0"),
                        x="x-axis",
                        y="frequency",
                        dodge=True,
                        color=pal,
                        alpha=0.8,
                        edgecolor="black",
                        linewidth=1,
                        order=x_axis_order,
                        ax=ax,
                        jitter=0.25,
                        size=4,
                        marker=marker,
                    )
                if ax.legend_:
                    ax.legend_.remove()
    control_df = manipulate_df(combined_oas_df, key=["oas_subject", "locus", "oas_author"])
    control_df = (
        control_df.query(f"distance_to_known_lcdr3 == {thresh}")
        .assign(dose_group="OAS 5 Len")
        .rename({"locus": "locus_light"}, axis=1)
    )

    for kappa_lambda in ["IGK", "IGL"]:
        if kappa_lambda == "IGK":
            marker = "o"
        else:
            marker = "^"
        sns.boxplot(
            data=control_df.query("frequency > 0").query(f"locus_light == '{kappa_lambda}'"),
            x="dose_group",
            y="frequency",
            hue="locus_light",
            dodge=True,
            palette={"IGK": "#F2F0F2", "IGL": "#F2F0F2"},
            linewidth=1,
            hue_order=["IGK", "IGL"],
            ax=distance_control,
            whis=[10, 90],
            fliersize=0,
        )
        adjust_boxplot(distance_control)
        sns.stripplot(
            data=control_df.query(f"locus_light == '{kappa_lambda}'"),
            x="dose_group",
            y="frequency",
            hue="locus_light",
            dodge=True,
            palette={"IGK": "#F2F0F2", "IGL": "#F2F0F2"},
            edgecolor="black",
            linewidth=0.5,
            alpha=0.8,
            size=4,
            hue_order=["IGK", "IGL"],
            ax=distance_control,
            marker=marker,
        )
        distance_control.legend_.remove()
    low_dose_label = r"20 $\mu$g"
    high_dose_label = r"100 $\mu$g"
    for index, (ax, title) in enumerate(
        [
            (distance_vrc01_combined, "VRC01 Class"),
            (distance_nonvrc01_combined, "Non-VRC01 Class\nwith 5aa LCDR3"),
            (distance_control, "Control"),
        ]
    ):
        if index in [0, 1]:
            x_ticks = [
                low_dose_label + "\nIGK",
                low_dose_label + "\nIGL",
                high_dose_label + "\nIGK",
                high_dose_label + "\nIGL",
            ]
            ax.set_xticklabels(
                x_ticks,
                rotation=0,
                fontsize=7,
            )
        else:
            ax.set_xticklabels(["OAS 5aa L3"], rotation=0, fontsize=7)
        if ax.legend_:
            ax.legend_.remove()
        ax.set_xlabel("")
        if index != 0:
            ax.set_yticklabels([])
            ax.set_ylabel("")
        else:
            label = "% BCRs with VRC01-class\nbnAb LCDR3"
            ax.set_ylabel(label, fontsize=8)
            ax.tick_params(axis="y", labelsize=8)

        ax.yaxis.set_major_locator(mtick.MultipleLocator(0.2))
        ax.yaxis.set_minor_locator(mtick.MultipleLocator(0.1))
        if index == 0:
            ax.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1, symbol=None))
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_ylim(-0.04, 1.04)
        if len(title.split("\n")) > 1:
            ax.set_title(title)
        else:
            ax.set_title(title, pad=10)

    distance_legend.set_xticks([])
    distance_legend.set_yticks([])
    distance_legend.axis("off")
    custom_lines = []
    custom_lines.append(
        Line2D([0], [0], markersize=4, markerfacecolor="white", color="black", linewidth=0, marker="o", label="IGK")
    )
    custom_lines.append(
        Line2D([0], [0], markersize=4, markerfacecolor="white", color="black", linewidth=0, marker="^", label="IGL")
    )
    distance_control.legend(
        handles=custom_lines,
        loc="upper left",
        frameon=False,
        handlelength=0.8,
        ncol=1,
        fontsize=8,
        bbox_to_anchor=(0.8, 1.05),
        labelspacing=0.1,
    )


def plot_key_mutations(
    fig, gs, palette: dict, mutational_df: pd.DataFrame, vrc01_ref_airr_df: pd.DataFrame, cottrell_path: Path
):
    """Have Key mutations"""
    panel_e_gs = gs
    key_low_dose = fig.add_subplot(panel_e_gs[0, 0])
    key_high_dose = fig.add_subplot(panel_e_gs[0, 1])
    key_control = fig.add_subplot(panel_e_gs[0, 2])

    # get mutational sets of control VH1-2 antibodies
    lat = LinkedAirrTable(vrc01_ref_airr_df.assign(is_vrc01_class=True))
    controal_lat_with_mutations = add_mutational_sets(lat, vrc01_ref_airr_df, cottrell_path)

    # get mutational set for our target antibodies
    mutational_df["sequence_id_heavy"] = mutational_df.index
    mutational_df["sequence_id_light"] = mutational_df.index
    lat = LinkedAirrTable(
        mutational_df[~((mutational_df["cellid"].isna()) | (mutational_df["mutations_heavy"].isna()))].reset_index(
            drop=True
        ),
        key_column="cellid",
    )
    lat_with_mutational_sets = add_mutational_sets(lat, vrc01_ref_airr_df, cottrell_path)
    plottable_df = lat_with_mutational_sets.query("weeks_post > 0").query("vaccine_group=='vaccine'").copy()

    # y = residues comes from groupby function
    y = "residues"
    key = [
        "pubid",
        "timepoint",
        "dose_group",
        "weeks_post",
        "vaccine_group",
        "is_vrc01_class",
    ]
    plottable_df_grouped = pd.DataFrame(plottable_df).groupby(key).apply(get_top_n_percet).reset_index()
    for df, ax, color in [
        (
            plottable_df_grouped.query("dose_group=='Low Dose'").query("is_vrc01_class"),
            key_low_dose,
            palette["Low Dose"]["vrc01_class_kappa"],
        ),
        (
            plottable_df_grouped.query("dose_group=='High Dose'").query("is_vrc01_class"),
            key_high_dose,
            palette["High Dose"]["vrc01_class_kappa"],
        ),
    ]:
        sns.boxplot(
            data=df,
            x="weeks_post",
            y=y,
            dodge=False,
            color=color,
            order=df.sort_values("weeks_post")["weeks_post"].unique(),
            ax=ax,
            linewidth=1,
            whis=[10, 90],
            fliersize=0,
        )
        adjust_boxplot(ax)
        sns.stripplot(
            data=df,
            x="weeks_post",
            y=y,
            color=color,
            dodge=False,
            linewidth=1,
            alpha=0.8,
            edgecolor="black",
            order=df.sort_values("weeks_post")["weeks_post"].unique(),
            size=4,
            ax=ax,
        )

    # control frame has no groupby since each value is a bNAb
    y = "cottrell_focused_v_common_score"
    sns.boxplot(
        data=controal_lat_with_mutations,
        x="is_vrc01_class",
        y=y,
        color="#F2F0F2",
        linewidth=1,
        ax=key_control,
        fliersize=0,
        whis=[10, 90],
    )
    adjust_boxplot(ax)
    sns.stripplot(
        data=controal_lat_with_mutations,
        x="is_vrc01_class",
        y=y,
        color="#F2F0F2",
        linewidth=1,
        alpha=0.8,
        edgecolor="black",
        size=4,
        ax=key_control,
    )
    # controal_lat_with_mutations.to_csv("controal_lat_with_mutations.csv")
    key_control.set_xticklabels(["VRC01 Class\nbnAbs"])
    low_dose_label = r"20 $\mu$g Dose"
    high_dose_label = r"100 $\mu$g Dose"

    for index, ax in enumerate([key_low_dose, key_high_dose, key_control]):
        if index in [0, 1]:
            if index == 0:
                ax.set_ylabel("Number of key VRC01-class\nresidues in heavy chain", fontsize=8, labelpad=5)
            if index == 0:
                ax.set_title(low_dose_label)
            else:
                ax.set_title(high_dose_label)
            ax.set_xlabel("Weeks Post Vaccination")
            new_labels = format_xtick_labels(ax.get_xticklabels())
            ax.set_xticklabels(new_labels, fontsize=7)
        elif index == 2:
            ax.set_xlabel("")
            ax.set_title("Control")
        if index == 1:
            ax.set_yticklabels([])
        if index > 0:
            ax.set_ylabel("")
        if index == 0 or index == 1:
            ax.set_ylim(0, 4.2)
        else:
            ax.set_ylim(0, 18)
        if index == 2:
            ax.yaxis.set_major_locator(mtick.MultipleLocator(4))
            ax.yaxis.set_minor_locator(mtick.MultipleLocator(2))
        else:
            ax.yaxis.set_major_locator(mtick.MultipleLocator(2))
            ax.yaxis.set_minor_locator(mtick.MultipleLocator(1))

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        if ax.legend_:
            ax.legend_.remove()


def plot_vh12_non_5_len_distro(data: Data, outpath: str):
    """
    Paper 2, LCDR3 length distribution of non VRC01 class VH12
    """
    unblind = data.get_unblinded_df().dropna(how="any", subset=["cdr3_aa_light"])
    unblind["cdr3_aa_light_length"] = unblind["cdr3_aa_light"].str.len().astype(int)
    both_vrc01_and_non = unblind.query("`has_vh1-2`==True").query("timepoint!='V02'").query("vaccine_group!='Placebo'")
    dekosky = data.get_dekosky()
    dekosky["cdr3_aa_light_length"] = dekosky["cdr3_aa_light"].str.len().astype(int)
    dekosky["locus_light"] = dekosky["v_call_top_light"].str[0:3]
    palette = data.get_plot_parameters().get_pallete()
    p = {"IGK": palette["Low Dose"]["vrc01_class_kappa"], "IGL": palette["High Dose"]["vrc01_class_lambda"]}
    print(palette)
    _, row_ax1 = plt.subplots(
        3,
        1,
        figsize=(4.5, 9),
        sharex=False,
        sharey=False,
    )

    flat_ax = row_ax1.flatten()
    for index, (df, ax) in enumerate(zip([both_vrc01_and_non, dekosky], flat_ax[:2])):
        groupby = "locus_light"
        plottable_df = df.groupby([groupby]).apply(lambda x: x["cdr3_aa_light_length"].value_counts(normalize=True))
        if isinstance(plottable_df, pd.Series):
            plottable_df = plottable_df.reset_index().rename(
                {"level_1": "cdr3_aa_light_length", "cdr3_aa_light_length": "count"}, axis=1
            )
        else:
            plottable_df = plottable_df.stack().reset_index().rename({0: "count"}, axis=1)
        sns.barplot(
            data=plottable_df,
            palette=p,
            x="cdr3_aa_light_length",
            y="count",
            ax=ax,
            hue=groupby,
            order=range(3, 14),
            linewidth=1,
            edgecolor="black",
        )
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_ylim(0, 1)
        if index == 1:
            ax.set_ylabel("Frequency", fontsize=12, labelpad=7)
            ax.legend_.set_title("Locus")
            ax.legend_.remove()
            ax.set_xlabel("")
            ax.set_title("DeKosky IGHV1-2", fontsize=14)
        if index == 0:
            ax.set_ylabel("Frequency", fontsize=12, labelpad=7)
            ax.legend(title="", loc="upper right", fontsize="small", fancybox=True)
            # ax.legend_.remove()
            ax.set_xlabel("")
            ax.set_title("G001 IGHV1-2", fontsize=14)

        ax.yaxis.set_minor_locator(mtick.MultipleLocator(0.1))
        # ax.set_xlabel("LCDR3 Length")

    #'enrichment
    lambda_dekosky = dekosky[dekosky["locus_light"] == "IGL"]["cdr3_aa_light_length"].value_counts(normalize=True)
    lambda_us = both_vrc01_and_non.query("locus_light=='IGL'")["cdr3_aa_light_length"].value_counts(normalize=True)
    df1 = (lambda_us / lambda_dekosky).reset_index().rename({"index": "length"}, axis=1).assign(locus="lambda")
    kappa_dekosky = dekosky[dekosky["locus_light"] == "IGK"]["cdr3_aa_light_length"].value_counts(normalize=True)
    kappa_us = both_vrc01_and_non.query("locus_light=='IGK'")["cdr3_aa_light_length"].value_counts(normalize=True)
    df2 = (kappa_us / kappa_dekosky).reset_index().rename({"index": "length"}, axis=1).assign(locus="kappa")
    enrich_df = pd.concat([df1, df2])
    p = {"kappa": palette["Low Dose"]["vrc01_class_kappa"], "lambda": palette["High Dose"]["vrc01_class_lambda"]}
    ax = sns.barplot(
        data=enrich_df,
        x="length",
        y="cdr3_aa_light_length",
        palette=p,
        hue="locus",
        edgecolor="black",
        hue_order=["kappa", "lambda"],
        ax=flat_ax[-1],
        order=range(3, 14),
    )
    ax.set_ylabel("Ratio", fontsize=12)
    ax.legend_.remove()
    ax.set_xlabel("LCDR3 Length", fontsize=12)
    ax.set(yscale="log")
    ax.set_title("Ratio G001/DeKosky", fontsize=14)
    sns.despine()
    # ax.set_xlim(-0.5, 8.5)
    plt.subplots_adjust(left=0.15, top=0.95, wspace=0.05, bottom=0.15, hspace=0.5)
    plt.savefig(outpath + ".svg")
    plt.savefig(outpath + ".pdf")
    plt.savefig(outpath + ".png", dpi=300)


def plot_temporal_sequence_features(data: Data, outpath: str):
    """
    Sequence features of VRc01 and non vrc01 but stratifed by time. S8 in paper 2
    """
    unblinded_df = data.get_unblinded_df()
    mutational_unblinded_df = data.get_mutational_df()
    oas_5_len = data.get_oas_5_len_df()
    vrc01_class_df = data.get_vrc01_ref_airr()
    oas_vh12_df = data.get_oas_vh12_df()
    plot_params = data.get_plot_parameters()

    plottable_df, _, _, oas_vh12_df = _add_sequence_metrics(
        unblinded_df, mutational_unblinded_df, oas_5_len, vrc01_class_df, oas_vh12_df
    )
    # use subplots, not gridspec
    fig, (row_ax1, row_ax2, row_ax3, row_ax4, row_ax5) = plt.subplots(
        5,
        5,
        figsize=(8.5, 11),
        gridspec_kw={"width_ratios": [1, 1, 1, 1, 0.2]},
        sharex=False,
        sharey=False,
    )

    # row 1 light chains
    low_dose_label = r"20 $\mu$g"
    high_dose_label = r"100 $\mu$g"
    titles_map = {"Low Dose": low_dose_label, "High Dose": high_dose_label}
    class_map = {False: "non-VRC01-class", True: "VRC01-class"}
    light_chain_df = plottable_df.groupby(["dose_group", "is_vrc01_class", "weeks_post"]).apply(
        lambda x: x["light_plot"].value_counts(normalize=True)
    )
    light_chain_df = light_chain_df.reset_index().rename({"level_3": "light_plot", "light_plot": "frequency"}, axis=1)
    for index, (ax, (is_vrc01, dose_group)) in enumerate(
        zip(row_ax1[:-1], product([True, False], ["Low Dose", "High Dose"]))
    ):

        df = light_chain_df.query(f"is_vrc01_class=={is_vrc01}").query(f"dose_group=='{dose_group}'")
        pivoted_plot = df.pivot("weeks_post", "light_plot", "frequency").fillna(0.0).sort_index()
        colors = [explicit_mapping[i] for i in sort_order]
        pivoted_plot.loc[:, sort_order].plot(
            kind="bar",
            stacked=True,
            color=colors,
            width=0.75,
            linewidth=1,
            edgecolor="black",
            ax=ax,
        )
        if index != 0:
            ax.set_yticklabels([])
            ax.set_ylabel("")
        else:
            label = "% BCRs with VRC01-class \n" + r"bnAb $\mathregular{V_{K/L}}$"
            ax.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1, symbol=None))
            ax.set_ylabel(label)
        order = sorted(plottable_df["weeks_post"].unique())
        ax.set_xticklabels(order, rotation=0)
        ax.yaxis.set_major_locator(mtick.MultipleLocator(0.2))
        ax.yaxis.set_minor_locator(mtick.MultipleLocator(0.1))
        ax.legend_.remove()
        ax.set_ylim(0, 1.05)
        ax.set_xlabel("")
        ax.set_title(f"{titles_map[dose_group]} {class_map[is_vrc01]}")
    sns.despine()
    custom_lines = []
    for label in sort_order:
        color = explicit_mapping[label]
        if label != "Other":
            crude_label = label[2:].replace("V", "")
        else:
            crude_label = label
        custom_lines.append(Patch(facecolor=color, edgecolor="black", linewidth=1, label=crude_label)),
    row_ax1[-1].legend(
        custom_lines,
        [i._label for i in custom_lines],
        loc="lower center",
        frameon=False,
        handlelength=0.8,
        ncol=1,
        fontsize=8,
        bbox_to_anchor=(0.5, -0.05),
        labelspacing=0.1,
    )
    row_ax1[-1].axis("off")

    # row 2, 5 Q and E and -5 W
    pal = plot_params.get_pallete()
    pal = [
        pal["Low Dose"]["vrc01_class"],
        pal["High Dose"]["vrc01_class"],
        pal["Low Dose"]["nonvrc01_class"],
        pal["High Dose"]["nonvrc01_class"],
    ]
    for sup_ax, metric in zip([row_ax2, row_ax5], ["has_QE_at_96", "has_w100b"]):

        def _determine_metric(df: pd.DataFrame) -> pd.Series:
            return pd.Series({True: df[metric].sum() / len(df), False: 1 - (df[metric].sum() / len(df))})

        if metric == "has_QE_at_96":
            sub_df = plottable_df.query("has_5_len_lcdr3")
        else:
            sub_df = plottable_df.query("`has_vh1-2`")

        combined = (
            sub_df.groupby(["dose_group", "pubid", "is_vrc01_class", "weeks_post"])
            .apply(_determine_metric)
            .stack()
            .reset_index()
            .rename({"level_4": metric, 0: "freq"}, axis=1)
        ).query(metric)
        order = sorted(combined["weeks_post"].unique())
        for index, (ax, (is_vrc01, dose_group)) in enumerate(
            zip(sup_ax[:-1], product([True, False], ["Low Dose", "High Dose"]))
        ):

            color = pal[index]
            df = combined.query(f"is_vrc01_class=={is_vrc01}").query(f"dose_group=='{dose_group}'")
            sns.boxplot(
                data=df,
                x="weeks_post",
                y="freq",
                linewidth=1,
                color=color,
                ax=ax,
                whis=[10, 90],
                fliersize=0,
                order=order,
            )
            sns.stripplot(
                data=df,
                x="weeks_post",
                y="freq",
                linewidth=1,
                color=color,
                ax=ax,
                jitter=0.15,
                size=4,
                order=order,
                edgecolor="black",
            )
            adjust_boxplot(ax)
            if index != 0:
                ax.set_yticklabels([])
                ax.set_ylabel("")
            else:
                ax.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1, symbol=None))
                if metric == "has_QE_at_96":
                    y_label = "% BCRs with E/Q at\nLCDR3 position 96"
                else:
                    y_label = "% BCRs with\n" + r"$\mathregular{Trp_{105-3}}$ in HCDR3"
                ax.set_ylabel(y_label)
            ax.yaxis.set_major_locator(mtick.MultipleLocator(0.2))
            ax.yaxis.set_minor_locator(mtick.MultipleLocator(0.1))
            ax.set_ylim(-0.04, 1.05)
            if metric == "has_w100b":
                ax.set_xlabel("Weeks Post Vaccine")
            else:
                ax.set_xlabel("")
        sup_ax[-1].axis("off")
    key = [
        "pubid",
        "dose_group",
        "locus_light",
        "is_vrc01_class",
        "vaccine_group",
        "weeks_post",
    ]
    class_per_timepoint = manipulate_df(plottable_df.query("has_5_len_lcdr3").copy(), key=key)
    order = sorted(plottable_df["weeks_post"].unique())
    for thresh, sup_ax in zip([0, 1], [row_ax3, row_ax4]):
        for index, (ax, (is_vrc01, dose_group)) in enumerate(
            zip(sup_ax[:-1], product([True, False], ["Low Dose", "High Dose"]))
        ):
            for kappa_lambda in ["IGK", "IGL"]:
                plot_df = class_per_timepoint.query(
                    f"locus_light == '{kappa_lambda}' and dose_group == '{dose_group}' and is_vrc01_class =={is_vrc01} and distance_to_known_lcdr3 == {thresh}"
                )
                color = pal[index]
                color = {"IGK": color, "IGL": color}
                if kappa_lambda == "IGK":
                    marker = "o"
                else:
                    marker = "^"
                if not plot_df.query("frequency > 0").empty:
                    sns.boxplot(
                        data=plot_df.query("frequency > 0"),
                        x="weeks_post",
                        y="frequency",
                        dodge=True,
                        linewidth=1,
                        hue="locus_light",
                        palette=color,
                        ax=ax,
                        whis=[10, 90],
                        order=order,
                        fliersize=0,
                        hue_order=["IGK", "IGL"],
                    )
                    sns.stripplot(
                        data=plot_df.query("frequency > 0"),
                        x="weeks_post",
                        y="frequency",
                        dodge=True,
                        edgecolor="black",
                        hue="locus_light",
                        linewidth=1,
                        palette=color,
                        order=order,
                        ax=ax,
                        size=4,
                        hue_order=["IGK", "IGL"],
                        marker=marker,
                    )
                    adjust_boxplot(ax)
                if not plot_df.query("frequency  == 0").empty:
                    sns.stripplot(
                        data=plot_df.query("frequency == 0"),
                        x="weeks_post",
                        y="frequency",
                        dodge=True,
                        hue="locus_light",
                        palette=color,
                        alpha=0.8,
                        order=order,
                        edgecolor="black",
                        linewidth=1,
                        ax=ax,
                        jitter=0.25,
                        size=4,
                        hue_order=["IGK", "IGL"],
                        marker=marker,
                    )
                if index != 0:
                    ax.set_yticklabels([])
                    ax.set_ylabel("")
                else:
                    ax.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1, symbol=None))
                    if thresh == 0:
                        y_label = "% BCRs with VRC01-class\nbnAb LCDR3"
                    else:
                        y_label = "% BCRs with VRC01-class\nbnAb LCDR3" + r"$\leq$" + str(thresh) + " mut"
                    ax.set_ylabel(y_label)
                ax.set_ylim(-0.04, 1.05)
                ax.set_xlabel("")
                # ax.set_xlabel("Weeks Post Vaccine")
                ax.yaxis.set_major_locator(mtick.MultipleLocator(0.2))
                ax.yaxis.set_minor_locator(mtick.MultipleLocator(0.1))
                if ax.legend_:
                    ax.legend_.remove()
        custom_lines = []
        custom_lines.append(
            Line2D([0], [0], markersize=4, markerfacecolor="white", color="black", linewidth=0, marker="o", label="IGK")
        )
        custom_lines.append(
            Line2D([0], [0], markersize=4, markerfacecolor="white", color="black", linewidth=0, marker="^", label="IGL")
        )
        sup_ax[-1].legend(
            handles=custom_lines,
            loc="upper left",
            frameon=False,
            handlelength=0.8,
            ncol=1,
            fontsize=8,
            bbox_to_anchor=(0.01, 1.05),
            labelspacing=0.1,
        )
        sup_ax[-1].axis("off")
    sns.despine()
    plt.subplots_adjust(hspace=0.17, wspace=0.05, top=0.98, right=0.97, left=0.08, bottom=0.05)
    fig.savefig(outpath + ".svg")
    fig.savefig(outpath + ".png", dpi=300)


def plot_seq_logos_by_gene(data: Data, outpath: str):
    """plot seq logos by specific v call, S18 of paper 2"""
    unblinded_df = data.get_unblinded_df()
    unblined_vrc01 = unblinded_df.query("is_vrc01_class == True")
    kappa_vrc01 = data.get_kappa_vrc01_select_df()
    lambda_vrc01 = data.get_lambda_vrc01_select_df()
    human_naive = data.get_human_naive_5_len_df().query("is_vrc01_class")
    fig, all_axis = plt.subplots(6, 2, figsize=(6, 4), sharex=False, sharey=False)

    # do kappa first
    # get all first column axis
    col_1 = list(map(lambda x: x[0], all_axis))
    cdr3_aa_light_list = kappa_vrc01["cdr3_aa_no_gaps"].to_list()
    _plot_sequence_logo(col_1[0], True, ytitle="bnAbs", char_set=cdr3_aa_light_list, xticks=[], title="IGK")
    for v_gene, ax in zip(["IGKV1-33", "IGKV3-20", "IGKV1-5", "IGKV3-15"], col_1[1:-1]):
        v_gene_df = unblined_vrc01[
            unblined_vrc01["v_call_light"].str.split(",").apply(lambda x: any([v_gene in i for i in x]))
        ]
        xticks = []
        cdr3_aa_light_list = v_gene_df["cdr3_aa_light"].to_list()
        _plot_sequence_logo(ax, True, ytitle=v_gene, char_set=cdr3_aa_light_list, xticks=xticks)

    cdr3_aa_light_list = human_naive.query("locus_light=='IGK'")["cdr3_aa_light"].to_list()
    xticks = list(range(93, 98))
    _plot_sequence_logo(col_1[-1], True, ytitle="Human naive\nGT8 binders", char_set=cdr3_aa_light_list, xticks=xticks)
    col_1[-1].set_xlabel("Light Chain Position")

    # now do lambda, just plot it one at a time wiht second column axis
    col_2 = list(map(lambda x: x[1], all_axis))
    cdr3_aa_light_list = lambda_vrc01["cdr3_aa_no_gaps"].to_list()

    # start plotting at first axis in row
    row_start = 0
    _plot_sequence_logo(col_2[row_start], True, ytitle="bnAbs", char_set=cdr3_aa_light_list, xticks=[], title="IGL")
    collect = []
    for c in ["IGLV2-14", "IGLV2-23", "IGLV2-11"]:
        v_gene_df = unblined_vrc01[
            unblined_vrc01["v_call_light"].str.split(",").apply(lambda x: any([c in i for i in x]))
        ]
        collect.append(v_gene_df)
    v_gene_df = pd.concat(collect)
    cdr3_aa_light_list = v_gene_df["cdr3_aa_light"].to_list()
    _plot_sequence_logo(col_2[row_start + 1], True, ytitle="Select IGL", char_set=cdr3_aa_light_list, xticks=[])

    cdr3_aa_light_list = human_naive.query("locus_light=='IGL'")["cdr3_aa_light"].to_list()
    xticks = list(range(93, 98))
    _plot_sequence_logo(
        col_2[row_start + 2], True, ytitle="Human naive\nGT8 binders", char_set=cdr3_aa_light_list, xticks=xticks
    )
    col_2[row_start + 2].set_xlabel("Light Chain Position")
    if row_start == 0:
        for ax in col_2[row_start + 3 :]:
            ax.axis("off")
    else:
        for ax in col_2[0 : row_start + 3]:
            ax.axis("off")

    plt.subplots_adjust(
        top=0.93,
        left=0.15,
        right=0.9,
        hspace=0.05,
        wspace=0.70,
        bottom=0.15,
    )
    fig.savefig(outpath + ".svg")
    fig.savefig(outpath + ".png", dpi=300)


def plot_heavy_chain_seq_logos_and_jchain(data: Data, outpath: str):
    unblinded_df = data.get_unblinded_df().query("Plate_Type_heavy == 'Probe Specific'")
    human_naive = data.get_human_naive_5_len_df()
    bnabs = data.get_vrc01_ref_airr().assign(is_vrc01_class=True)
    fig, all_axis = plt.subplots(
        4, 4, figsize=(8.5, 5), gridspec_kw={"width_ratios": [1, 0.6, 1, 0.6]}, sharex=False, sharey=False
    )
    j_genes = ["IGHJ1", "IGHJ2", "IGHJ3", "IGHJ4", "IGHJ5", "IGHJ6"]
    low_dose_label = r"20 $\mu$g"
    high_dose_label = r"100 $\mu$g"

    # unpack first column as a list of axis
    x_ticks = list(range(-5, 0, 1))
    low_dose = unblinded_df.query("dose_group == 'Low Dose'")
    high_dose = unblinded_df.query("dose_group == 'High Dose'")
    for col_index in range(0, 4, 2):
        # vrc01 class
        col_1 = list(map(lambda x: x[col_index], all_axis))

        # vrc01 class j gene distro
        col_2 = list(map(lambda x: x[col_index + 1], all_axis))
        for row_inx, (title, dataframe) in enumerate(
            [
                ("bnAbs", bnabs),
                (low_dose_label, low_dose),
                (high_dose_label, high_dose),
                ("Human\nnaive\nGT8\nbinders", human_naive),
            ]
        ):
            if row_inx > 0:
                if col_index == 0:
                    is_vrc01_class = True
                else:
                    is_vrc01_class = False
            else:
                is_vrc01_class = True
            plottable_df = dataframe.query(f"is_vrc01_class=={is_vrc01_class}")
            cdr3_list = plottable_df["cdr3_aa_heavy"].str[-5:].to_list()
            # if col_index > 1:
            #     title = None
            if col_index == 2 and row_inx == 0:
                title = "VRC01\nclass\nbnAbs"
            _plot_sequence_logo(
                col_1[row_inx],
                True,
                ytitle=title,
                char_set=cdr3_list,
                xticks=x_ticks,
                ytitle_font=10,
            )
            if row_inx == 3:
                col_1[row_inx].set_xlabel("HCDR3 Position")
            (
                plottable_df["j_call_heavy"]
                .str.split(",")
                .str.get(0)
                .str.split("*")
                .str.get(0)
                .value_counts(normalize=True)
                .reindex(j_genes)
                .fillna(0)
                .sort_index()
                * 100
            ).plot(ax=col_2[row_inx], kind="bar", color="grey", ylim=(0, 60))
            col_2[row_inx].set_ylabel("%")
            if row_inx == 3:
                col_2[row_inx].set_xlabel("IGHJ Gene")
            if row_inx == 0 and col_index == 0:
                col_1[row_inx].annotate("VRC01 Class", xy=(1, 1.2), ha="center", xycoords="axes fraction")
            elif row_inx == 0 and col_index == 2:
                col_1[row_inx].annotate("Non-VRC01 Class", xy=(1, 1.2), ha="center", xycoords="axes fraction")
    plt.subplots_adjust(
        top=0.93,
        left=0.07,
        right=0.95,
        hspace=0.13,
        wspace=0.3,
        bottom=0.15,
    )
    sns.despine()
    fig.savefig(outpath + ".svg")
    fig.savefig(outpath + ".png", dpi=300)


def plot_indel_features(data: Data, outpath: str):
    mutational_df = data.get_mutational_df()
    palette = data.plot_parameters.get_pallete()

    low_dose_label = r"20 $\mu$g"
    high_dose_label = r"100 $\mu$g"

    def get_lcdr1_feature(X, feature):
        features_in_lcdr1 = 0
        # print(X)
        if X is None:
            return 0

        for x in X:
            position = int(re.findall(r"\d+", x)[0])
            if position > 26 and position < 33:
                if x[-1] == feature:
                    features_in_lcdr1 += 1
        return features_in_lcdr1

    def _calculate_indels(df: pd.DataFrame) -> pd.DataFrame:
        "Calculate the number of indel and for heavy and light chain"
        df["heavy_del"] = df["v_sequence_alignment_aa_heavy"].str.count("-")
        df["light_del"] = df["v_sequence_alignment_aa_light"].str.count("-")
        df["heavy_ins"] = df["v_germline_alignment_aa_heavy"].str.count("-")
        df["light_ins"] = df["v_germline_alignment_aa_light"].str.count("-")
        df["num_del_lcdr1"] = df["mutations_light"].apply(lambda x: get_lcdr1_feature(x, "-"))
        df["num_r_lcdr3"] = df["mutations_light"].apply(lambda x: get_lcdr1_feature(x, "R"))
        return df

    def _calculate_indel_bool(df: pd.DataFrame) -> pd.DataFrame:
        """Calculate if there is an indel in the dataframe of interest (i.e. a subset dataframe"""
        df["has_heavy_del"] = df["heavy_del"] > 0
        df["has_light_del"] = df["light_del"] > 0
        df["has_heavy_ins"] = df["heavy_ins"] > 0
        df["has_light_ins"] = df["light_ins"] > 0
        df["has_lcdr1_del"] = df["num_del_lcdr1"] > 0
        df["has_lcdr1_R"] = df["num_r_lcdr3"] > 0
        return df

    # get a dataframe with indel counts and a bool if there is an indel at all
    post_vaccine_with_indels = _calculate_indel_bool(
        _calculate_indels(mutational_df)
        .query("Plate_Type_heavy=='Probe Specific'")
        .query("timepoint!='V02'")
        .query("vaccine_group!='Placebo'")
    )

    # group by first and second shot
    post_vaccine_with_indels.loc[
        post_vaccine_with_indels[post_vaccine_with_indels["weeks_post"] <= 8].index,
        "shot",
    ] = "first shot"
    post_vaccine_with_indels.loc[
        post_vaccine_with_indels[post_vaccine_with_indels["weeks_post"] > 8].index, "shot"
    ] = "second shot"

    fig, axes = plt.subplots(6, 2, figsize=(7, 7))
    row_labels = [
        # r"%$\mathregular{V_{K1-33}}$" + "LCDR1\nmutations to R",
        # r"%$\mathregular{V_{K/L}}$" + "LCDR1\nmutations to R",
        r"%$\mathregular{V_{K3-20}}$ LCDR1" + "\ndeletions",
        r"%$\mathregular{V_{K3-20}}$" + "\ndeletions",
        r"%$\mathregular{V_{K/L}}$" + "\ndeletions",
        r"%$\mathregular{V_{H}}$" + "\ndeletions",
        r"%$\mathregular{V_{K/L}}$" + "\ninsertions",
        r"%$\mathregular{V_{H}}$" + "\ninsertions",
    ]
    row_ranges = [
        # (-0.01, 0.3),
        # (-0.01, 0.3),
        (-0.01, 0.17),
        (-0.01, 0.17),
        (-0.01, 0.11),
        (-0.01, 0.17),
        (-0.001, 0.05),
        (-0.01, 0.12),
    ]
    for index, (row_ax, row_label, row_range, metric) in enumerate(
        zip(
            axes,
            row_labels,
            row_ranges,
            [
                # "has_lcdr1_R",
                # "has_lcdr1_R",
                "has_lcdr1_del",
                "has_light_del",
                "has_light_del",
                "has_heavy_del",
                "has_light_ins",
                "has_heavy_ins",
            ],
        )
    ):

        # if index == 0:
        #     plotting_df = post_vaccine_with_indels[
        #         post_vaccine_with_indels["v_call_top_light"].str.split("*").str.get(0) == "IGKV1-33"
        #     ].copy()
        if index < 2:
            plotting_df = post_vaccine_with_indels[
                post_vaccine_with_indels["v_call_top_light"].str.split("*").str.get(0) == "IGKV3-20"
            ].copy()

        else:
            plotting_df = post_vaccine_with_indels.copy()

        def _determine_metric(df: pd.DataFrame) -> pd.Series:
            return pd.Series({True: df[metric].sum() / len(df), False: 1 - (df[metric].sum() / len(df))})

        combined = (
            plotting_df.groupby(["dose_group", "pubid", "is_vrc01_class", "shot"])
            .apply(_determine_metric)
            .stack()
            .reset_index()
            .rename({"level_4": metric, 0: "freq"}, axis=1)
        )
        combined["x_axis"] = combined["dose_group"] + "_" + combined["shot"]
        combined_vrc01 = combined.query("is_vrc01_class == True").query(metric)
        combined_nonvrc01 = combined.query("is_vrc01_class != True").query(metric)
        for col_ax, vrc01_class in zip(row_ax, [True, False]):
            x_nest = ["Low Dose_first shot", "Low Dose_second shot", "High Dose_first shot", "High Dose_second shot"]
            for shot in ["first shot", "second shot"]:
                for dose in ["Low Dose", "High Dose"]:
                    x_axis = dose + "_" + shot
                    if vrc01_class:
                        c = "vrc01_class"
                        df = combined_vrc01
                    else:
                        c = "nonvrc01_class"
                        df = combined_nonvrc01
                    color = palette[dose][c]
                    sns.boxplot(
                        data=df.query(f"x_axis == '{x_axis}'"),
                        x="x_axis",
                        y="freq",
                        linewidth=1,
                        color=color,
                        ax=col_ax,
                        whis=[10, 90],
                        order=x_nest,
                        fliersize=0,
                    )
                    sns.stripplot(
                        data=df.query(f"x_axis == '{x_axis}'"),
                        x="x_axis",
                        y="freq",
                        linewidth=1,
                        color=color,
                        ax=col_ax,
                        order=x_nest,
                        size=4,
                        edgecolor="black",
                    )
            if index == 4:
                col_ax.yaxis.set_major_locator(mtick.MultipleLocator(0.01))
            else:
                col_ax.yaxis.set_major_locator(mtick.MultipleLocator(0.05))
            col_ax.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1, symbol=None, decimals=0))
            # col_ax.yaxis.set_major_formatter(mtick.FuncFormatter(lambda _, x: 100 * x))
            adjust_boxplot(col_ax)
            col_ax.set_xticklabels([])
            col_ax.set_xlabel("")
            col_ax.set_ylabel("")
            col_ax.set_ylim(row_range)
        if index == 4:
            row_ax[0].set_ylabel(row_label, va="bottom", labelpad=10)
        else:
            row_ax[0].set_ylabel(row_label, va="bottom")

    x_ticks = [
        low_dose_label + "\nPrime",
        low_dose_label + "\nBoost",
        high_dose_label + "\nPrime",
        high_dose_label + "\nBoost",
    ]
    axes[-1, 0].set_xticklabels(x_ticks)
    axes[-1, 1].set_xticklabels(x_ticks)
    axes[0, 0].set_title("VRC01-class")
    axes[0, 1].set_title("non-VRC01-class")

    sns.despine()
    plt.tight_layout()
    fig.savefig(outpath + ".svg")
    fig.savefig(outpath + ".png", dpi=300)


def plot_temporal_w_check(data: Data, outpath: str):
    """plot temporal trp check at -5 position , S10 of paper 2

    Parameters
    ----------
    data : Data
        [description]
    outpath : str
        [description]

    Returns
    -------
    [type]
        [description]
    """
    # check W in HCDR3 as a function of 5 length
    unblinded_df = data.get_unblinded_df()
    mutational_unblinded_df = data.get_mutational_df()
    oas_5_len = data.get_oas_5_len_df()
    vrc01_class_df = data.get_vrc01_ref_airr()
    oas_vh12_df = data.get_oas_vh12_df()
    plot_params = data.get_plot_parameters()
    pallete = plot_params.get_pallete()
    plottable_df, _, _, oas_vh12_df = _add_sequence_metrics(
        unblinded_df, mutational_unblinded_df, oas_5_len, vrc01_class_df, oas_vh12_df
    )
    # use subplots, not gridspec
    fig, (row_ax1) = plt.subplots(1, 2, figsize=(7, 3), sharex=True, sharey=True)

    # row 1 light chains
    metric = "has_w100b"
    low_dose_label = r"20 $\mu$g Dose"
    high_dose_label = r"100 $\mu$g Dose"
    titles_map = {"Low Dose": low_dose_label, "High Dose": high_dose_label}
    sub_df = plottable_df.query("has_5_len_lcdr3").query("is_vrc01_class == False")

    def _determine_metric(df: pd.DataFrame) -> pd.Series:
        return pd.Series({True: df[metric].sum() / len(df), False: 1 - (df[metric].sum() / len(df))})

    combined = (
        sub_df.groupby(["dose_group", "pubid", "weeks_post"])
        .apply(_determine_metric)
        .stack()
        .reset_index()
        .rename({"level_3": "has_w100b", 0: "freq"}, axis=1)
    ).query("has_w100b")
    pal = [pallete["Low Dose"]["nonvrc01_class"], pallete["High Dose"]["nonvrc01_class"]]
    order = sorted(plottable_df["weeks_post"].unique())
    for index, (ax, dose) in enumerate(zip([row_ax1[0], row_ax1[1]], ["Low Dose", "High Dose"])):
        df = combined.query(f"dose_group == '{dose}'")
        color = pal[index]
        sns.boxplot(
            data=df,
            x="weeks_post",
            y="freq",
            linewidth=1,
            color=color,
            ax=ax,
            whis=[10, 90],
            fliersize=0,
            order=order,
        )
        sns.stripplot(
            data=df,
            x="weeks_post",
            y="freq",
            linewidth=1,
            color=color,
            ax=ax,
            jitter=0.15,
            size=4,
            order=order,
            edgecolor="black",
        )
        adjust_boxplot(ax)
        ax.set_ylim(-0.04, 1.05)
        if index == 0:
            ax.set_ylabel("Frequency")
            ax.yaxis.set_major_locator(mtick.MultipleLocator(0.2))
            ax.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1, symbol=None))
            ax.yaxis.set_minor_locator(mtick.MultipleLocator(0.1))
        else:
            ax.set_ylabel("")
        ax.set_xlabel("Weeks Post Vaccine")
        ax.set_title(titles_map[dose])
    label = "% BCRs with " + r"$\mathregular{Trp_{105-3}}$" + " in HCDR3 among non-VRC01-class BCRs with 5-aa LCDR3s"
    fig.suptitle(label)
    sns.despine()
    plt.tight_layout()
    fig.savefig(outpath + ".svg")
    fig.savefig(outpath + ".png", dpi=300)


def plot_combined_metrics(data: Data, outpath: str):
    # combined metrics
    unblinded_df = data.get_unblinded_df()
    mutational_unblinded_df = data.get_mutational_df()
    oas_5_len = data.get_oas_5_len_df()
    vrc01_class_df = data.get_vrc01_ref_airr()
    oas_vh12_df = data.get_oas_vh12_df()
    plottable_df, _, _, oas_vh12_df = _add_sequence_metrics(
        unblinded_df, mutational_unblinded_df, oas_5_len, vrc01_class_df, oas_vh12_df
    )
    plottable_df["distance_to_known_lcdr3"] = plottable_df[
        ["distance_to_known_lcdr3_kappa", "distance_to_known_lcdr3_lambda"]
    ].min(axis=1)
    plottable_df["has_distance"] = plottable_df["distance_to_known_lcdr3"] <= 1
    plottable_df["has_light_segment"] = plottable_df["light_plot"] != "Other"
    plottable_df["features"] = plottable_df[["has_distance", "has_light_segment", "has_w100b", "has_QE_at_96"]].sum(
        axis=1
    )
    vrc01_only_plot = plottable_df.query("is_vrc01_class == True")

    def _determine_metric(df):
        return df["features"].value_counts(normalize=True).reindex(range(0, 5)).fillna(0).sort_index()[::-1].cumsum()

    sub_plot = (
        vrc01_only_plot.groupby(["dose_group", "pubid", "shot"])
        .apply(_determine_metric)
        .stack()
        .reset_index()
        .rename({0: "freq"}, axis=1)
    )
    sub_plot["x-axis"] = sub_plot["dose_group"] + sub_plot["shot"]
    _, ax = plt.subplots(1, 1, figsize=(7, 3.0))
    order = ["Low Dosefirst shot", "Low Dosesecond shot", "High Dosefirst shot", "High Dosesecond shot"]
    col = ["#6C65FF", "#FF65AB", "#F8FF65", "#65FFB9"]
    pal = dict(zip(range(1, 5), col))
    sns.boxplot(
        data=sub_plot,
        x="x-axis",
        y="freq",
        hue="features",
        linewidth=1,
        ax=ax,
        palette=pal,
        hue_order=[1, 2, 3, 4],
        whis=[10, 90],
        fliersize=0,
        order=order,
    )
    sns.stripplot(
        data=sub_plot,
        x="x-axis",
        y="freq",
        hue="features",
        dodge=True,
        linewidth=1,
        palette=pal,
        ax=ax,
        size=5,
        hue_order=[1, 2, 3, 4],
        order=order,
        edgecolor="black",
    )
    ax.legend_.remove()
    low_dose_label = r"20 $\mu$g"
    high_dose_label = r"100 $\mu$g"
    ax.set_xticklabels(
        [
            low_dose_label + "\nPrime",
            low_dose_label + "\nBoost",
            high_dose_label + "\nPrime",
            high_dose_label + "\nBoost",
        ],
        rotation=0,
        fontsize=10,
    )
    adjust_boxplot(ax)
    ax.set_xlabel("")
    ax.set_ylabel("Frequency of 1, 2, 3, or 4\nbnAb features in VRC01-class BCRs", fontsize=10)
    ax.tick_params(labelsize=12)
    ax.yaxis.set_major_locator(mtick.MultipleLocator(0.2))
    ax.yaxis.set_minor_locator(mtick.MultipleLocator(0.1))
    ax.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1))
    custom_lines = []
    for features, column in zip([1, 2, 3, 4], col):
        custom_lines.append(
            Line2D(
                [0],
                [0],
                markersize=4,
                markerfacecolor=column,
                color="black",
                linewidth=0,
                marker="o",
                label=f"{features}",
            )
        )
    legend = ax.legend(
        handles=custom_lines,
        loc="center left",
        frameon=False,
        handlelength=0.8,
        ncol=1,
        fontsize=10,
        bbox_to_anchor=(1.0, 0.5),
        labelspacing=0.1,
        title="Number of\nFeatures",
    )

    plt.setp(legend.get_title(), multialignment="center")
    sns.despine()
    plt.tight_layout()
    plt.savefig(outpath + ".pdf")
    plt.savefig(outpath + ".svg")
    plt.savefig(outpath + ".png", dpi=300)


def plot_multi_week_comparison_vh_mut(data: Data, outpath: str, week_one: int, week_two: int):
    unblinded_df = data.get_unblinded_df()
    plot_params = data.get_plot_parameters()
    full_palette = plot_params.get_pallete()

    # only get vaccinated and FNA
    unblinded_df_vac = (
        unblinded_df.query("vaccine_group=='vaccine'")
        .query(f"weeks_post == {week_one} or weeks_post == {week_two}")
        .reset_index(drop=True)
    )
    unblinded_df_vac["v_mutation_heavy"] = unblinded_df_vac["v_mutation_heavy"] * 100
    unblinded_df_vac["v_mutation_light"] = unblinded_df_vac["v_mutation_light"] * 100

    def get_donors_with(is_vrc01):
        """Only get donors that have both timepoints"""
        ubdonor = (
            unblinded_df_vac.query(f"is_vrc01_class=={is_vrc01}")
            .groupby("PubID")
            .apply(lambda x: len(x["timepoint"].unique()) > 1)
            .reset_index()
        )
        donors_with_multi = ubdonor[ubdonor[0]]["PubID"].to_list()
        return donors_with_multi

    plottable = (
        unblinded_df_vac[unblinded_df_vac["PubID"].isin(get_donors_with(True))]
        .reset_index(drop=True)
        .query("is_vrc01_class")
    )
    total_donors = len(plottable["PubID"].unique())
    num_rows = math.ceil(total_donors / 6)
    fig, axes = plt.subplots(num_rows, 6, figsize=(8.5, 2.0 + (1.5 * num_rows - 1)), sharex=False, sharey=True)
    gb_keys = sorted(
        plottable.groupby(["dose_group", "PubID"]).groups.keys(),
        key=lambda x: x[0],
        reverse=True,
    )
    gb_df = plottable.groupby(["dose_group", "PubID"])
    last_index = 0
    for gb_key, (ax_index, ax) in zip(gb_keys, enumerate(axes.flatten())):
        g_df = gb_df.get_group(gb_key)
        dg = g_df["dose_group"].iloc[0]
        c = full_palette[dg]["vrc01_class"]
        ptid = gb_key[1].split("_")[1]
        sns.boxplot(x="weeks_post", y="v_mutation_heavy", ax=ax, data=g_df, whis=[10.0, 90.0], fliersize=0, color=c)
        sns.stripplot(
            x="weeks_post", y="v_mutation_heavy", ax=ax, data=g_df, size=5, edgecolor="black", linewidth=1, color=c
        )
        ax.set_title(ptid)
        # ax.set_yticklabels([])
        if ax_index == 0:
            label = r"$\mathregular{V_H}$" + " gene\n%mutation (nt)"
            ax.set_ylabel(label)
        else:
            ax.set_ylabel("")
        adjust_boxplot(ax)
        ax.set_xlabel("")
        last_index = ax_index
    for ax in axes.flatten()[last_index + 1 :]:
        ax.axis("off")
    sns.despine()
    plt.tight_layout()
    fig.savefig(outpath + "A.svg")
    fig.savefig(outpath + "A.png", dpi=300)

    # B will be non-VRC01
    plottable = (
        unblinded_df_vac[unblinded_df_vac["PubID"].isin(get_donors_with(False))]
        .reset_index(drop=True)
        .query("is_vrc01_class==False")
    )
    total_donors = len(plottable["PubID"].unique())
    num_rows = math.ceil(total_donors / 6)
    fig, axes = plt.subplots(num_rows, 6, figsize=(8.5, 2.0 + (1.5 * num_rows - 1)), sharex=False, sharey=True)
    axes = axes.flatten()
    gb_keys = sorted(
        plottable.groupby(["dose_group", "PubID"]).groups.keys(),
        key=lambda x: x[0],
        reverse=True,
    )
    gb_df = plottable.groupby(["dose_group", "PubID"])
    last_index = 0
    for gb_key, (ax_index, ax) in zip(gb_keys, enumerate(axes)):
        g_df = gb_df.get_group(gb_key)
        dg = g_df["dose_group"].iloc[0]
        c = full_palette[dg]["nonvrc01_class"]
        ptid = gb_key[1].split("_")[1]
        sns.boxplot(x="weeks_post", y="v_mutation_heavy", ax=ax, data=g_df, whis=[10, 90], fliersize=0, color=c)
        sns.stripplot(
            x="weeks_post", y="v_mutation_heavy", ax=ax, data=g_df, size=5, edgecolor="black", linewidth=1, color=c
        )
        ax.set_title(ptid)
        if ax_index % 6 == 0:
            label = r"$\mathregular{V_H}$" + " gene\n%mutation (nt)"
            ax.set_ylabel(label)
        else:
            ax.set_ylabel("")
        adjust_boxplot(ax)
        ax.set_xlabel("")
        last_index = ax_index
    for ax in axes[last_index + 1 :]:
        ax.axis("off")
    sns.despine()
    plt.tight_layout()
    fig.savefig(outpath + "B.svg")
    fig.savefig(outpath + "B.png", dpi=300)
