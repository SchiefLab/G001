from g001.data import Data
from g001.figures.util import MinorSymLogLocator, adjust_boxplot
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch
import pandas as pd
import seaborn as sns
from g001.figures.somatic import format_xtick_labels


def plot_BCR_vs_cluster(
    cluster_df: pd.DataFrame,
    palette: dict,
    is_vrc01_class: bool,
    row_ax: list[plt.Axes],
    xlabel: str = "Number of VRC01-class BCR Seqs",
):
    row_ax_1, row_ax_2 = row_ax[0], row_ax[1]
    plottable_df = cluster_df.groupby(["pubid", "dose_group", "vaccine_group"], as_index=False).apply(
        lambda x: pd.Series({"num_clusters": len(x["cluster"].unique()), "num_seqs": len(x)})
    )
    for index, (dose, ax) in enumerate(zip(palette, [row_ax_1, row_ax_2])):
        color = palette[dose]
        sns.scatterplot(
            data=plottable_df.query(f"dose_group=='{dose}'"),
            x="num_seqs",
            y="num_clusters",
            edgecolor="black",
            color=color,
            s=40,
            linewidth=1,
            ax=ax,
        )
        if is_vrc01_class:
            ax.set_xlim(-10, 260)
            ax.set_ylim(-10, 260)
            X = range(-10, 300)
        else:
            ax.set_xlim(-10, 400)
            ax.set_ylim(-10, 400)
            X = range(-10, 400)
        if not is_vrc01_class and dose == "High Dose":
            color = "gold"
        Y = X
        ax.plot(
            X,
            Y,
            linestyle="--",
            color=color,
            zorder=0,
            markeredgecolor="black",
            # edgecolor="black",
            # z_order=0,
            # alpha=0.5,
        )
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        if index == 0:
            ax.set_title("Low Dose")
        else:
            # ax.set_yticklabels([])
            # ax.set_ylabel("")
            ax.set_title("High Dose")
        ax.set_ylabel("Number of Clusters")
        ax.set_xlabel(xlabel)


def summarize_cluster(cluster_df: pd.DataFrame, is_vrc01_class: bool) -> pd.DataFrame:
    antibodies = cluster_df.query(f"is_vrc01_class=={is_vrc01_class}")
    results = (
        pd.DataFrame(antibodies)
        .groupby(["cluster"])
        .apply(
            lambda x: pd.Series(
                {
                    "size": len(x),
                    "subjects": list(x["pubid"].unique()),
                    "num_subjects": len(list(x["pubid"].unique())),
                    "timepoints": list(x["timepoint"].unique()),
                    "num_timeponts": len(list(x["timepoint"].unique())),
                    "vdj_mut_aa_mean": x["v_mutation_aa_heavy"].mean(),
                    "vdj_mut_aa_std": x["v_mutation_aa_heavy"].std(),
                    "cdr3_aa_heavy": x["cdr3_aa_heavy"].to_list(),
                }
            )
        )
        .sort_values("num_subjects")[::-1]
    )
    return results


def plot_summary_cluster(cluster_df: pd.DataFrame, is_vrc01_class: bool, palette: dict, ax: plt.Axes):
    cluster_df = cluster_df.query(f"is_vrc01_class=={is_vrc01_class}")
    summary_df = summarize_cluster(cluster_df, is_vrc01_class).reset_index()
    df_info = cluster_df.groupby(["dose_group", "cluster"]).head(1)[["cluster", "dose_group", "is_vrc01_class"]]
    cluster_summary_df = df_info.merge(summary_df, on="cluster")
    plottable_df = []
    for dose in ["Low Dose", "High Dose"]:
        dose_specific = cluster_summary_df.query(f"dose_group=='{dose}'")
        number_singletons = len(dose_specific.query("num_timeponts==1 and num_subjects==1")["cluster"].unique())
        number_timepoints = len(dose_specific.query("num_timeponts > 1")["cluster"].unique())
        number_donors = len(dose_specific.query("num_subjects > 1")["cluster"].unique())

        plottable_df += [
            {
                "count": number_singletons,
                "Dose_Group": dose,
                "label": "cluster_single_timepoint",
            },
            {
                "count": number_timepoints,
                "label": "cluster_more_one_timepoint",
                "Dose_Group": dose,
            },
            {
                "count": number_donors,
                "Dose_Group": dose,
                "label": "cluster_multiple_subjects",
            },
        ]
    plottable_df = pd.DataFrame(plottable_df)
    sns.barplot(
        data=plottable_df,
        hue="Dose_Group",
        y="count",
        x="label",
        # orient="h",
        edgecolor="black",
        ax=ax,
        palette=palette,
    )

    if is_vrc01_class:
        ax.set_ylim(0, 1100)
    else:
        pass
        ax.set_ylim(0, 10000)
    ax.set_yscale("symlog", linthresh=10)
    ax.set_xticklabels(
        ["Single Donor\n and Timepoint", "Multiple\nTimepoints", "Multiple\nDonors"],
        rotation=90,
        ha="right",
        va="center",
        rotation_mode="anchor",
    )
    ax.set_ylabel("Number of Clusters", labelpad=0)
    ax.set_xlabel("")
    ax.legend_.remove()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.yaxis.set_minor_locator(MinorSymLogLocator(10))  # custom class to handle sym log ticks.
    if is_vrc01_class:
        ax.set_yticklabels([r"$0$", r"$10$", r"$10^2$", r"$10^3$"], va="top")
    else:
        ax.set_yticklabels([r"$0$", r"$10$", r"$10^2$", r"$10^3$", r"$10^4$"], va="top")

    custom_lines = []
    for label, color in zip(["Low Dose", "High Dose"], [palette["Low Dose"], palette["High Dose"]]):
        custom_lines.append(Patch(facecolor=color, edgecolor="black", linewidth=1, label=label)),
    ax.legend(
        custom_lines,
        [i._label for i in custom_lines],
        loc="upper center",
        frameon=False,
        handlelength=0.8,
        ncol=1,
        fontsize=8,
        bbox_to_anchor=(0.8, 1.05),
        labelspacing=0.25,
    )


def plot_distance_histogram(distance_df: pd.DataFrame, ax: plt.Axes, xmax: int):
    pall_list = ["#949494", "#56B4E9", "#DF8F05", "#D55E00"]
    palette = dict(zip(["all_v_all", "inter_centroid", "intra_cluster"], pall_list))
    maximum = xmax
    for method in distance_df.method.unique():
        sns.histplot(
            data=distance_df.query("method==@method").query(f"distance < {maximum + 1}"),
            x="distance",
            hue="method",
            bins=30,
            stat="density",
            palette=palette,
            discrete=True,
            hue_order=distance_df.method.unique(),
            multiple="dodge",
            ax=ax,
        )
    ax.set_xlim(0, maximum)
    ax.legend_.remove()
    ax.set_ylim(0, 0.2)
    ax.set_ylabel("Frequency")
    ax.set_xlabel("Distance between BCRs")
    ax.xaxis.set_major_locator(plt.MultipleLocator(5))
    ax.xaxis.set_minor_locator(plt.MultipleLocator(1))

    custom_lines = []
    for label, color in zip(
        ["All vs All", "Intra-Cluster", "Inter-Centroid"], [pall_list[0], pall_list[2], pall_list[1]]
    ):
        custom_lines.append(Patch(facecolor=color, edgecolor="black", linewidth=1, label=label)),
    ax.legend(
        custom_lines,
        [i._label for i in custom_lines],
        loc="upper center",
        frameon=False,
        handlelength=0.8,
        ncol=2,
        fontsize=10,
        bbox_to_anchor=(0.5, -0.3),
        labelspacing=0.25,
    )


def plot_timecourse_clonality(
    cluster_df: pd.DataFrame, palette: dict[str, str], is_vrc01_class: bool, row_ax: list[plt.Axes]
):
    row_ax_1, row_ax_2 = row_ax[0], row_ax[1]
    plottable_df = (
        cluster_df.groupby(["pubid", "timepoint", "weeks_post", "is_vrc01_class", "dose_group"])
        .apply(
            lambda x: pd.Series(
                {
                    "clonality": len(x["cluster"].unique()) / len(x),
                    "num_seqs": len(x),
                    "num_clusters": len(x["cluster"].unique()),
                }
            )
        )
        .reset_index()
    )
    for index, (dose, ax) in enumerate(zip(palette, [row_ax_1, row_ax_2])):
        color = palette[dose]
        sns.boxplot(
            data=plottable_df.query("dose_group==@dose").query("is_vrc01_class==@is_vrc01_class"),
            x="weeks_post",
            y="clonality",
            color=color,
            ax=ax,
            fliersize=0,
            linewidth=1,
            whis=[10, 90],
        )
        adjust_boxplot(ax)
        sns.stripplot(
            data=plottable_df.query("dose_group==@dose").query("is_vrc01_class==@is_vrc01_class"),
            x="weeks_post",
            y="clonality",
            color=color,
            ax=ax,
            linewidth=1,
            edgecolor="black",
            size=5,
            jitter=0.1,
        )
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_ylim(0, 1.05)
        if index == 0:
            ax.set_title("Low Dose")
        else:
            # ax.set_yticklabels([])
            # ax.set_ylabel("")
            ax.set_title("High Dose")
        ax.set_ylabel("Fraction Unique Clones")
        ax.set_xlabel("Weeks Post Vaccination")
        ax.set_xticklabels(format_xtick_labels(ax.get_xticklabels()))


def plot_seqs_vs_cluster(cluster_df: pd.DataFrame, is_vrc01_class: bool, ax: plt.Axes, offset: int = 100):
    cluster_df = cluster_df.query(f"is_vrc01_class=={is_vrc01_class}").sort_values("cluster_len").reset_index(drop=True)
    cluster_df.index = pd.Int64Index(range(1, len(cluster_df) + 1))
    X = list(cluster_df.index)
    Y = []
    for cluster_len, cluster_len_df in cluster_df.groupby("cluster_len"):
        # horrible  hack!
        for _, _ in cluster_len_df.groupby("cluster"):
            if not Y:
                last_value = 0
            else:
                last_value = Y[-1]
            Y += [last_value + 1] * int(cluster_len)
            # print(sub_cluster_len, Y)

    X_1, Y_1 = (
        X[0 : cluster_df.query("cluster_len==1").index[-1] - 1],
        Y[0 : cluster_df.query("cluster_len==1").index[-1] - 1],
    )

    X_2, Y_2 = (
        X[cluster_df.query("cluster_len==1").index[-1] - 1 : cluster_df.query("cluster_len==2").index[-1] - 1],
        Y[cluster_df.query("cluster_len==1").index[-1] - 1 : cluster_df.query("cluster_len==2").index[-1] - 1],
    )
    X_3, Y_3 = (
        X[cluster_df.query("cluster_len==2").index[-1] - 1 : cluster_df.query("cluster_len==3").index[-1] - 1],
        Y[cluster_df.query("cluster_len==2").index[-1] - 1 : cluster_df.query("cluster_len==3").index[-1] - 1],
    )
    X_4, Y_4 = (
        X[cluster_df.query("cluster_len==3").index[-1] - 1 : cluster_df.query("cluster_len>3").index[-1] - 1],
        Y[cluster_df.query("cluster_len==3").index[-1] - 1 : cluster_df.query("cluster_len>3").index[-1] - 1],
    )
    pall_list = ["#949494", "#56B4E9", "#DF8F05", "#D55E00"]
    for color, x, y in zip(pall_list, [X_1, X_2, X_3, X_4], [Y_1, Y_2, Y_3, Y_4]):
        ax.plot(x, y, color=color, linewidth=5)
    ax.hlines(X[-1], xmin=-100, xmax=X[-1], linestyles="--", color="grey")
    ax.hlines(Y_4[-1], xmin=-100, xmax=X[-1], linestyles="--", color="black")
    if is_vrc01_class:
        ax.set_xlim(-50, 3000)
        ax.set_ylim(-50, 3000)
    else:
        pass
    ax.hlines(X[-1], xmin=-100, xmax=X[-1], linestyles="--", color="grey")
    ax.hlines(Y_4[-1], xmin=-100, xmax=X[-1], linestyles="--", color="black")
    ax.set_ylabel("Number of Clusters")
    ax.set_xlabel("Number of BCR Sequences")
    sns.despine()
    custom_lines = []
    custom_lines.append(Patch(facecolor="white", edgecolor="white", linewidth=0.1, label="BCRs per cluster:"))
    custom_lines.append(Patch(facecolor="white", edgecolor="white", linewidth=1, label=""))
    for label, color in zip(["N=1", "N=2", "N=3", "N>3"], pall_list):
        custom_lines.append(Patch(facecolor=color, edgecolor="black", linewidth=1, label=label)),
    ax.legend(
        custom_lines,
        [i._label for i in custom_lines],
        loc="upper center",
        frameon=False,
        handlelength=0.8,
        ncol=3,
        fontsize=10,
        bbox_to_anchor=(0.35, -0.3),
        labelspacing=0.25,
    )
    ax.annotate(
        f"Total Clusters = {round((Y[-1]/X[-1])*100,1)}%\nof BCRs",
        xy=(X[-1] / 2, Y_4[-1] + offset),
        ha="center",
    )
    if is_vrc01_class:
        ax.annotate(f"VRC01-class BCRs={X[-1]}", xy=(X[-1] / 2, X[-1] + 100), ha="center")
    else:
        ax.annotate(f"nonVRC01-class BCRs={X[-1]}", xy=(X[-1] / 2, X[-1] + 100), ha="center")


def plot_cluster_size_histogram(
    cluster_df: pd.DataFrame, palette: dict[str, str], is_vrc01_class: bool, row_ax: list[plt.Axes]
):
    row_ax_1, row_ax_2 = row_ax[0], row_ax[1]
    plottable_df = (
        pd.DataFrame(cluster_df)
        .groupby(["dose_group", "cluster"], as_index=False)
        .apply(lambda x: pd.Series({"num_clusters": len(x)}))
    )
    for index, (dose, ax) in enumerate(zip(palette, [row_ax_1, row_ax_2])):
        color = palette[dose]
        sns.histplot(
            data=plottable_df.query(f"dose_group=='{dose}'"),
            x="num_clusters",
            edgecolor="black",
            linewidth=1,
            color=color,
            binwidth=1,
            ax=ax,
        )
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        if index == 0:
            ax.set_title("Low Dose")
        else:
            ax.set_title("High Dose")
        ax.set_xlabel("Cluster Size")
        ax.set_yscale("log")
        if is_vrc01_class:
            ax.set_ylim(0, 2000)
            ax.set_xlim(1, 80)
        else:
            # pass
            ax.set_ylim(0, 3000)
            ax.set_xlim(1, 80)

        ax.set_ylabel("Count")


def plot_class_clonality(data: Data, is_vrc01_class: bool = True) -> plt.figure:
    fig = plt.figure(figsize=(8.5, 11))
    cluster_df = data.get_unblinded_sequences().query(f"is_vrc01_class=={is_vrc01_class}").query("timepoint!='V02'")
    cluster_df = cluster_df[~cluster_df["cluster"].isna()].reset_index()
    distance_df = data.get_distance_df(is_vrc01_class)
    plot_parms = data.plot_parameters
    gs = GridSpec(13, 14, figure=fig, height_ratios=[1, 1, 1, 0.1, 1, 1, 1, 1, 1, 1, 1, 1, 1])

    palette = plot_parms.get_pallete()
    if is_vrc01_class:
        hue_tille = "vrc01_class"
    else:
        hue_tille = "nonvrc01_class"

    parsed_pal = {"Low Dose": palette["Low Dose"][hue_tille], "High Dose": palette["High Dose"][hue_tille]}

    # Row 1 Ax
    row1_ax_1 = fig.add_subplot(gs[0:3, :5])
    row1_ax_2 = fig.add_subplot(gs[0:3, 5:11])
    row1_ax_3 = fig.add_subplot(gs[0:3, 11:])
    if is_vrc01_class:
        offset = 100
    else:
        offset = 100
    plot_seqs_vs_cluster(cluster_df, is_vrc01_class, row1_ax_1, offset=offset)
    if is_vrc01_class:
        xmax = 30
    else:
        xmax = 55
    plot_distance_histogram(distance_df, row1_ax_2, xmax)
    plot_summary_cluster(cluster_df, is_vrc01_class, parsed_pal, row1_ax_3)

    # Row 2 Ax, cluster histogrram
    row2_ax_1 = fig.add_subplot(gs[4:7, :7])
    row2_ax_2 = fig.add_subplot(gs[4:7, 7:])
    plot_cluster_size_histogram(cluster_df, parsed_pal, is_vrc01_class, [row2_ax_1, row2_ax_2])

    # Row 3 Ax, BCR vs cluster
    row3_ax_1 = fig.add_subplot(gs[7:10, :7])
    row3_ax_2 = fig.add_subplot(gs[7:10, 7:])
    plot_BCR_vs_cluster(cluster_df, parsed_pal, is_vrc01_class, [row3_ax_1, row3_ax_2])

    # Row 4 Ax, timecourse clonality
    row4_ax_1 = fig.add_subplot(gs[10:13, :7])
    row4_ax_2 = fig.add_subplot(gs[10:13, 7:])
    plot_timecourse_clonality(cluster_df, parsed_pal, is_vrc01_class, [row4_ax_1, row4_ax_2])

    fig.subplots_adjust(top=0.94, wspace=20, hspace=50, right=0.96, left=0.09)
    return fig
