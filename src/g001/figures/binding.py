import matplotlib.ticker as mtick
import pandas as pd
import seaborn as sns
from matplotlib import gridspec
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from g001.data import Data
from g001.figures.tables import plot_binding_table
from g001.figures.util import MinorSymLogLocator, adjust_boxplot

log_fixed_locator_lookup = {
    -14: "10 fM",
    -13: "100 fM",
    -12: "1 pM",
    -11: "10 pM",
    -10: "100 pM",
    -9: "1 nM",
    -8: "10 nM",
    -7: "100 nM",
    -6: "1 $\mathregular{\mu}$M",  # noqa: W605
    -5: "10 $\mathregular{\mu}$M",  # noqa: W605
    -4: "$\mathregular{\geq}$100$\mathregular{\mu}$M",  # noqa: W605
    -3: "1 mM",
    -2: "10 mM",
    -1: "100 mM",
    0: "1 M",
    1: "10 M",
}


fixed_locator_lookup = {
    1e-13: "100 fM",
    1e-12: "1 pM",
    1e-11: "10 pM",
    1e-10: "100 pM",
    1e-09: "1 nM",
    1e-08: "10 nM",
    1e-07: "100 nM",
    1e-06: "1 $\mathregular{\mu}$M",  # noqa: W605
    1e-05: "10 $\mathregular{\mu}$M",  # noqa: W605
    0.0001: "$\mathregular{\geq}$100$\mathregular{\mu}$M",  # noqa: W605
    0.001: "1 mM",
    0.01: "10 mM",
    0.1: "100 mM",
    0: "1 M",
    1: "10 M",
    2: "100 M",
    10: "1 G",
}


def tens_formatter(x):
    if int(x) == 0:
        return "1"
    return "$10^{{{0:d}}}$".format(int(x))
    # return r"$\mathregular{10^{" + str(x) + r"}}$"


def tinv(p, df):
    # used for 95% confidence interval
    return abs(t.ppf(p / 2, df))


def get_total_mutations(X):
    heavy = (X["v_mutation_aa_heavy"]) * len(X["v_sequence_alignment_aa_heavy"])
    light = (X["v_mutation_aa_light"]) * len(X["v_sequence_alignment_aa_light"])
    total = round(heavy) + round(light)
    return int(total)


def plot_strip_on_box(
    ax,
    dataframe,
    palette,
    x="x-axis",
    y="KD_fix",
    hue="is_vrc01_class",
    sort_by="weeks_post",
    hue_order=[True, False],
    lt_criteria=True,
) -> None:
    # plot boxplot then strip with KDs < 100 uM
    sns.boxplot(
        x=x,
        y=y,
        data=dataframe,
        hue=hue,
        dodge=True,
        palette=palette,
        linewidth=1,
        order=dataframe.sort_values(sort_by)[x].unique(),
        hue_order=hue_order,
        ax=ax,
        whis=[10, 90],
        # width=0.95,
        fliersize=0,
    )
    if lt_criteria:
        data = dataframe.query(f"{y} <1e-4")
    else:
        data = dataframe

    sns.stripplot(
        x=x,
        y=y,
        data=data,
        hue=hue,
        dodge=True,
        edgecolor="black",
        linewidth=1,
        alpha=0.8,
        order=dataframe.sort_values(sort_by)[x].unique(),
        size=5,
        hue_order=hue_order,
        palette=palette,
        ax=ax,
    )
    sns.stripplot(
        x=x,
        y=y,
        data=dataframe.query("KD_fix ==1e-4"),
        hue=hue,
        dodge=True,
        edgecolor="black",
        linewidth=1,
        alpha=0.8,
        order=dataframe.sort_values(sort_by)[x].unique(),
        size=5,
        jitter=0.25,
        hue_order=hue_order,
        palette=palette,
        ax=ax,
    )


def plot_gt8_binding(data: Data):
    """Plot GT8 Binding for G001 Ligands"""
    mature_spr = data.get_mature_spr_df()
    clk_spr = data.get_clk_spr_df()
    igl_spr = data.get_igl_spr_df()
    spr_df = pd.concat([mature_spr, clk_spr, igl_spr]).reset_index(drop=True)
    plot_parm = data.plot_parameters
    palette = plot_parm.get_pallete()
    memory_timepoint = plot_parm.get_memory_visit_lookup()
    pb_timepoint = plot_parm.get_pb_visit_lookup()
    gc_timepoint = plot_parm.get_gc_visit_lookup()
    annotate_args = plot_parm.get_small_annotate_args()

    fig = plt.figure(figsize=(9, 5.5))
    title_xy = (0.5, 1.3)
    gs = gridspec.GridSpec(
        3,
        10,
        figure=fig,
        height_ratios=[0.15, 0.8, 0.10],
        hspace=0.05,
        wspace=0.25,
        left=0.13,
        top=0.93,
        right=0.98,
        bottom=0.05,
    )

    # table specs
    clk_table_spec = fig.add_subplot(gs[0, 0])
    ig_table_spec = fig.add_subplot(gs[0, 1])
    pbmc_table_spec = fig.add_subplot(gs[0, 2:7])
    gc_table_spec = fig.add_subplot(gs[0, 7:9])
    pb_table_spec = fig.add_subplot(gs[0, 9:])

    # plot A spec
    clk_spec = fig.add_subplot(gs[1, 0])
    ig_spec = fig.add_subplot(gs[1, 1])
    pbmc_spec = fig.add_subplot(gs[1, 2:7])
    gc_spec = fig.add_subplot(gs[1, 7:9])
    pb_spec = fig.add_subplot(gs[1, 9:])

    # label axis
    label_axis = fig.add_subplot(gs[2, :])

    # first only get GT8 analytes
    single_analyte_spr = spr_df.query("publish_analyte=='GT8'").reset_index(drop=True)

    # only have to change igl to say weeks_post = 0 and Timepoint == 0
    single_analyte_spr.loc[single_analyte_spr["is_igl"], "weeks_post"] = 0
    single_analyte_spr.loc[single_analyte_spr["is_igl"], "Timepoint"] = "V00"
    single_analyte_spr.loc[single_analyte_spr["is_igl"], "Specimen_Type"] = "iGL"
    single_analyte_spr = single_analyte_spr.astype({"weeks_post": int})
    single_analyte_spr["x-axis"] = single_analyte_spr[["weeks_post", "Specimen_Type"]].apply(
        lambda x: "week " + str(x[0]) + "\n" + x[1] if x[1] != "CLK" else "week 0", axis=1
    )

    # change the FNA to GC and PMBC to MBC
    single_analyte_spr["x-axis"] = single_analyte_spr["x-axis"].apply(
        lambda x: x.replace("FNA", "GC").replace("PBMC", "MBC") if x else x
    )

    box_and_strip = {True: palette["Low Dose"]["vrc01_class"], False: palette["Low Dose"]["nonvrc01_class"]}

    # plot main SPR plots for each analyte
    for ax_index, (table_ax, ax) in enumerate(
        zip(
            [clk_table_spec, ig_table_spec, pbmc_table_spec, gc_table_spec, pb_table_spec],
            [clk_spec, ig_spec, pbmc_spec, gc_spec, pb_spec],
        )
    ):
        timepoints = {}
        table_title = "Null"
        plot_table_row_title = False
        overide_stats = None
        if ax_index == 0:
            timepoints = {"CLK": 0}
            table_title = "Naive\nPrecursors"
            plot_table_row_title = True
            overide_stats = 12
        elif ax_index == 1:
            timepoints = {"V00": 0}
            table_title = "Inferred\nPrecursors"
        elif ax_index == 2:
            timepoints = memory_timepoint
            table_title = "Memory B Cells in PBMCs"
        elif ax_index == 3:
            timepoints = gc_timepoint
            table_title = "GC B Cells in LNs"
        elif ax_index == 4:
            timepoints = pb_timepoint
            table_title = "Plasmablasts\nin PBMCs"
        # Only get timepoints of intreset
        print(timepoints)  # to beat the linter
        sub_plot = single_analyte_spr.query("Timepoint in @timepoints.keys()").copy()
        if ax_index == 0:
            sub_plot["x-axis"] = "week 0"

        # plot table in catch all function
        plot_binding_table(
            table_ax, sub_plot, annotate_args, table_title, title_xy, plot_table_row_title, overide_stats
        )

        # plot strip plot on boxplot
        plot_strip_on_box(ax, sub_plot, box_and_strip)

        # remove legend, will add manually
        ax.legend_.remove()
        ax.set_ylim(1e-11, 0.00015)
        ax.set_xlabel("")
        ax.set(yscale="log")

        if ax_index == 0:
            ax.set_ylabel("$\mathregular{K_D}$", size=14)  # noqa: W605
            ax.yaxis.set_major_locator(mtick.LogLocator(base=10, numticks=9))
            # custom locator, shouldnt be linear anywhere
            ax.yaxis.set_minor_locator(MinorSymLogLocator(linthresh=1))
            ax.yaxis.set_major_formatter(mtick.FuncFormatter(lambda x, p: fixed_locator_lookup[x]))
        else:
            ax.set_ylabel("")
            ax.yaxis.set_major_locator(mtick.LogLocator(base=10, numticks=9))
            ax.yaxis.set_minor_locator(mtick.LogLocator(base=10, subs=range(0, 10), numticks=9))
            ax.yaxis.set_major_formatter(mtick.NullFormatter())
        adjust_boxplot(ax)
        ax.tick_params(which="major", length=5, labelsize=12, width=0.7)
        ax.tick_params(which="minor", length=3, width=0.5)
        ax.axhline(y=0.0001, linestyle="--", linewidth=1, color="black", zorder=-10)
    sns.despine()

    custom_lines = [
        Line2D([0], [0], color=palette["Low Dose"]["vrc01_class"], lw=6, label="VRC01"),
        Line2D([0], [0], color=palette["Low Dose"]["nonvrc01_class"], lw=6, label="Non-VRC01"),
    ]
    label_axis.legend(
        custom_lines,
        ["VRC01-class", "non-VRC01-class"],
        loc="lower center",
        frameon=False,
        ncol=3,
        fontsize=12,
        bbox_to_anchor=(0.5, -0.85),
    )
    label_axis.set_axis_off()
    return fig
