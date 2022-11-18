import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import matplotlib.ticker as mtick
from g001.data import Data


def shape_data(dataframe: pd.DataFrame) -> pd.DataFrame:
    """Find out responders over 1 shot, 2 shots and all shots stratified by treatment group"""

    # drop pre vax
    dataframe = dataframe.drop(dataframe[dataframe["Visit"] == "V02"].index)

    # get first shot and second shot by timepoint
    dataframe.loc[dataframe[dataframe["Visit"].isin(["V05", "V06", "V07"])].index, "shot"] = "first"
    dataframe.loc[dataframe[dataframe["Visit"].isin(["V07A", "V08", "V09", "V10"])].index, "shot"] = "two"

    # make a new dataframe so we can easily groupby
    new_df = []
    for g, g_df in dataframe.groupby(["Treatment", "shot"]):

        # responders is Response =1
        responders = len(g_df.query("Response==1")["PubID"].unique())

        # total is how many PubIDS are in this group
        total = len(g_df["PubID"].unique())

        # make a new dataframe with the groupby info
        new_df.append(
            {
                "Treatment": g[0],
                "Shot": g[1],
                "responders": responders,
                "total": total,
                "percentage_responders": responders / total,
            }
        )

    # have to do this once again so we can make a dataframe with both 1 and 2 shots
    for g, g_df in dataframe.groupby(["Treatment"]):
        responders = len(g_df.query("Response==1")["PubID"].unique())
        total = len(g_df["PubID"].unique())
        new_df.append(
            {
                "Treatment": g,
                "Shot": "all",
                "responders": responders,
                "total": total,
                "percentage_responders": responders / total,
            }
        )
    plotting_df = pd.DataFrame(new_df).sort_values(["Treatment", "Shot"])[::-1].reset_index(drop=True)
    return plotting_df


def plot_overview(data: Data):
    frequency_df = data.get_flow_and_frequency_data()
    plot_params = data.plot_parameters
    treatment_pallette = plot_params.get_treatment_pallete()
    left_title = "After 1st vaccination\n(Weeks 3, 4, 8)"
    middle_title = "After 2nd vaccination\n(Weeks 9, 10, 11, 16)"
    right_title = "After 1st or 2nd vaccination\n(All timepoints)"
    x_axis_rotation = 0

    # make a reshaped dataframe
    plotting_df = shape_data(frequency_df)

    def change_width(ax, new_value):
        for patch in ax.patches:
            current_width = patch.get_width()
            diff = current_width - new_value

            # we change the bar width
            patch.set_width(new_value)

            # we recenter the bar
            patch.set_x(patch.get_x() + diff * 0.5)

    # don't care about second shot for now, only first and all ahosts
    # plotting_df = plotting_df[plotting_df["Shot"] != "two"].copy()

    # basic col catplot with our custom palette
    g = sns.catplot(
        data=plotting_df,
        kind="bar",
        x="Treatment",
        y="percentage_responders",
        col="Shot",
        hue="Treatment",
        palette=treatment_pallette,
        dodge=False,
        edgecolor="black",
        height=3,
        aspect=0.85,
        linewidth=1.2,
    )

    # flatten axes to go through them
    for ax_index, ax in enumerate(g.axes.flat):

        change_width(ax, 0.75)
        # annotate each bar with the number of responders/totla
        for x in ax.get_xticklabels():

            # only need x position
            (x1, _) = x.get_position()

            # title is set by seaborn is Shot =
            title = ax.get_title().split("= ")[-1]
            lookup = x.get_text()

            # get the df of just the treatment we are looking at
            lookup_df = plotting_df[(plotting_df["Treatment"] == lookup) & (plotting_df["Shot"] == title)]
            if len(lookup_df) > 1:
                raise AssertionError(f"{lookup_df} is more than on")
            else:
                lookup_df = lookup_df.iloc[0]

            reponsers = lookup_df["responders"]
            total = lookup_df["total"]
            percentage = lookup_df["percentage_responders"]

            # Se the text with a 0.01 nudge
            ax.text(x1, percentage + 0.01, f"{reponsers}/{total}", fontsize=12, ha="center", fontweight="bold")

        # # Set the x labels for all axis
        # ax.set_xticklabels(
        #     ["Placebo", r"20 $\mu$g", r"100 $\mu$g"],
        #     rotation=x_axis_rotation,
        #     fontsize=14,
        # )

        if ax_index == 0:
            ax.yaxis.set_major_locator(mtick.MultipleLocator(base=0.1))
            ax.yaxis.set_major_formatter(mtick.FuncFormatter(lambda y, _: str(int(float(y) * 100))))
            # this is the only way I have found to set label size with formatter
            ax.tick_params(axis="y", labelsize=14)
            ax.set_title(left_title, fontweight="bold", pad=20, fontsize=12)
        if ax_index == 1:
            ax.set_title(middle_title, fontweight="bold", pad=20, fontsize=12)
        if ax_index == 2:
            ax.set_title(right_title, fontweight="bold", pad=20, fontsize=12)

        ax.set_xlabel("")
    g.set_ylabels("% VRC01-class responders", fontsize=14)
    g.set(ylim=(0, 1.0))
    g.set_xticklabels(
        ["Placebo", r"20 $\mu$g", r"100 $\mu$g"],
        rotation=x_axis_rotation,
        fontsize=14,
    )
    plt.tight_layout()
    return g.figure
