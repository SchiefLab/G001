import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from numpy import concatenate


def format_xtick_labels(L):
    "Format the xtick labels to have the week and sample type"
    new_labels = []
    lookup = {
        "-4": "wk -4\nMBC",
        "3": "wk 3\nGC",
        "4": "wk 4\nMBC",
        "8": "wk 8\nMBC",
        "9": "wk 9\nPB",
        "10": "wk 10\nMBC",
        "11": "wk 11\nGC",
        "16": "wk 16\nMBC",
    }
    for x in L:
        text = x.get_text()
        new_text = lookup[text]
        new_labels.append(new_text)
    return new_labels


def shape_df(df: pd.DataFrame) -> pd.DataFrame:
    """Get median mutation rates per donor per timepoint"""
    df = df.groupby(
        ["pubid", "timepoint", "vaccine_group", "dose_group", "is_vrc01_class", "weeks_post"], as_index=False
    )[["v_mutation_aa_heavy", "v_mutation_aa_light", "v_mutation_heavy", "v_mutation_light"]].median()
    df = df.drop(df[df["vaccine_group"] != "vaccine"].index)
    return df


def plot_somatic_mutation_frequencies(
    dataframe: pd.DataFrame,
    palette: dict[str, str],
    vrc01: str,
    molecule: str,
    sup_title: str,
    strip_size: int = 7,
) -> plt.figure:
    """Plot somatic mutation frequencies

    Parameters
    ----------
    dataframe : pd.DataFrame
        unblinded df with just the BCR of interest
    palette : dict
        color palette for the plot
    vrc01 : str
        vrc01 class, vrc01 or nonvrc01
    molecule : str
        nt or aa
    sup_title : str
        Title to place above figure
    strip_size : int, optional
        size of markers, by default 7

    Raises
    ------
    ValueError
        If molecule is not nt or aa
    """

    # Get the dataframe with median mutation rates per donor per timepoint
    plotting_df = shape_df(dataframe.query('vaccine_group=="vaccine"'))
    plotting_df["v_mutation_aa_heavy"] = plotting_df["v_mutation_aa_heavy"] * 100
    plotting_df["v_mutation_aa_light"] = plotting_df["v_mutation_aa_light"] * 100
    if molecule == "aa" and vrc01 == "nonvrc01":
        plotting_df.loc[plotting_df[plotting_df["v_mutation_aa_heavy"] >= 22].index, "v_mutation_aa_heavy"] = 22

    if molecule == "nt":
        plotting_df["v_mutation_heavy"] = plotting_df["v_mutation_heavy"] * 100
        plotting_df["v_mutation_light"] = plotting_df["v_mutation_light"] * 100
    # setup subplots for a 2x2 grid, figsize could become a variable
    fig, (row_ax1, row_ax2) = plt.subplots(
        2,
        2,
        figsize=(8.5, 5.5),
        gridspec_kw={
            "width_ratios": [1, 1],
            "height_ratios": [1, 1],
        },
        sharex=True,
        sharey=True,
    )  # type: ignore

    # if we want vrc01 class, parse that df
    if vrc01 == "vrc01":
        _class = plotting_df[plotting_df["is_vrc01_class"]].copy()
    elif vrc01 == "nonvrc01":
        _class = plotting_df[~plotting_df["is_vrc01_class"]].copy()
    else:
        ValueError(f"vrc01 must be either 'vrc01' or 'nonvrc01' but is {vrc01}")

    # divide low and high dose into two sep dataframes
    low_dose = _class.query("dose_group == 'Low Dose'")
    high_dose = _class.query("dose_group == 'High Dose'")

    # depending on what molceule is asked for is what row we look up
    if molecule == "nt":
        start = "v_mutation"
    elif molecule == "aa":
        start = "v_mutation_aa"
    else:
        raise ValueError("need aa or nt")

    # iterate trhough low,high, heavy and light chain per figure
    for ax_index, (df, y_label, ax) in enumerate(
        [
            (low_dose, f"{start}_heavy", row_ax1[0]),
            (high_dose, f"{start}_heavy", row_ax1[1]),
            (low_dose, f"{start}_light", row_ax2[0]),
            (high_dose, f"{start}_light", row_ax2[1]),
        ]
    ):
        sns.boxplot(
            ax=ax,
            data=df,
            x="weeks_post",
            y=y_label,
            fliersize=0,
            whis=False,
            hue="dose_group",
            palette=palette,
        )

        # adjust the boxplots to highlight median
        sns.stripplot(
            ax=ax,
            data=df,
            x="weeks_post",
            y=y_label,
            edgecolor="black",
            linewidth=1,
            palette=palette,
            hue="dose_group",
            marker="o",
            s=strip_size,
        )

        # y labels depending on which row we are on, VH or VK/VL with subscripts and moeldcule
        if ax_index == 0:
            label = r"$\mathregular{V_H}$" + f" gene\n%mutation ({molecule})"
            ax.set_ylabel(label)
        elif ax_index == 2:
            label = r"$\mathregular{V_{K/L}}$" + f" gene\n%mutation ({molecule})"
            ax.set_ylabel(label)

        # second column turns off labels
        else:
            ax.set_ylabel("")

        # for some reason this has to come first or it messes everything up
        ax.set_xlabel("")

        # get new labels with wk and sample type
        new_labels = format_xtick_labels(ax.get_xticklabels())
        ax.set_xticklabels(new_labels)

        # add the tick locator:
        if molecule == "aa":
            if vrc01 == "vrc01":
                ax.set_ylim(-0.5, 12.5)
            else:
                ax.set_ylim(-0.5, 23)
            ax.yaxis.set_major_locator(plt.MultipleLocator(2))
            ax.yaxis.set_minor_locator(plt.MultipleLocator(1))
            ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: str(int(y)) if y != 22 else r"$\geq22$"))
        else:
            if vrc01 == "vrc01":
                ax.set_ylim(0, 6.5)
            else:
                ax.set_ylim(0, 23)
            ax.yaxis.set_major_locator(plt.MultipleLocator(2))
            ax.yaxis.set_minor_locator(plt.MultipleLocator(1))

        # set column tittles that are dosage groups
        if ax_index == 0:
            ax.set_title(r"20 $\mu$g Dose", pad=20, fontsize=12)
        if ax_index == 1:
            ax.set_title(r"100 $\mu$g Dose", pad=20, fontsize=12)

    # finally remove spnes, not sure if this is necesary, but it is necessary to remove all legends
    for _, ax in enumerate(concatenate((row_ax1, row_ax2))):
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.get_legend().remove()

    # add a supfigure title slightly off centered
    fig.suptitle(sup_title, ha="center", x=0.53)
    fig.tight_layout()
    return fig


def plot_somatic_mutation_frequencies_violin(
    dataframe: pd.DataFrame,
    palette: dict[str, str],
    vrc01: str,
    molecule: str,
    sup_title: str,
    strip_size: int = 7,
) -> plt.figure:
    """Plot somatic mutation frequencies

    Parameters
    ----------
    dataframe : pd.DataFrame
        unblinded df with just the BCR of interest
    palette : dict
        color palette for the plot
    vrc01 : str
        vrc01 class, vrc01 or nonvrc01
    molecule : str
        nt or aa
    sup_title : str
        Title to place above figure
    outpath : str
        path to save figure
    strip_size : int, optional
        size of markers, by default 7

    Raises
    ------
    ValueError
        If molecule is not nt or aa
    """

    # Get the dataframe with median mutation rates per donor per timepoint
    plotting_df = dataframe.copy()
    order = sorted(list(plotting_df["weeks_post"].unique()))
    plotting_df["v_mutation_aa_heavy"] = plotting_df["v_mutation_aa_heavy"] * 100
    plotting_df["v_mutation_aa_light"] = plotting_df["v_mutation_aa_light"] * 100
    plotting_df["v_mutation_heavy"] = plotting_df["v_mutation_heavy"] * 100
    plotting_df["v_mutation_light"] = plotting_df["v_mutation_light"] * 100

    # setup subplots for a 2x2 grid, figsize could become a variable
    fig, (row_ax1, row_ax2) = plt.subplots(
        2,
        2,
        figsize=(8.5, 5.5),
        gridspec_kw={
            "width_ratios": [1, 1],
            "height_ratios": [1, 1],
        },
        sharex=True,
        sharey=True,
    )

    # if we want vrc01 class, parse that df
    if vrc01 == "vrc01":
        _class = plotting_df[plotting_df["is_vrc01_class"]].query("weeks_post != -4").copy()
    elif vrc01 == "nonvrc01":
        _class = plotting_df[~plotting_df["is_vrc01_class"]].copy()
    else:
        ValueError(f"vrc01 must be either 'vrc01' or 'nonvrc01' but is {vrc01}")

    # divide low and high dose into two sep dataframes
    low_dose = _class.query("dose_group == 'Low Dose'")
    high_dose = _class.query("dose_group == 'High Dose'")

    # depending on what molceule is asked for is what row we look up
    if molecule == "nt":
        start = "v_mutation"
    elif molecule == "aa":
        start = "v_mutation_aa"
    else:
        raise ValueError("need aa or nt")

    # iterate trhough low,high, heavy and light chain per figure
    for ax_index, (df, y_label, ax) in enumerate(
        [
            (low_dose, f"{start}_heavy", row_ax1[0]),
            (high_dose, f"{start}_heavy", row_ax1[1]),
            (low_dose, f"{start}_light", row_ax2[0]),
            (high_dose, f"{start}_light", row_ax2[1]),
        ]
    ):
        sns.violinplot(
            ax=ax,
            data=df,
            x="weeks_post",
            y=y_label,
            hue="dose_group",
            palette=palette,
            inner="quartile",
            saturation=1,
            scale="area",
            order=order,
            cut=0,
        )

        # y labels depending on which row we are on, VH or VK/VL with subscripts and moeldcule
        if ax_index == 0:
            label = r"$\mathregular{V_H}$" + f" gene\n%mutation ({molecule})"
            ax.set_ylabel(label)
        elif ax_index == 2:
            label = r"$\mathregular{V_{K/L}}$" + f" gene\n%mutation ({molecule})"
            ax.set_ylabel(label)

        # second column turns off labels
        else:
            ax.set_ylabel("")

        # for some reason this has to come first or it messes everything up
        ax.set_xlabel("")

        # get new labels with wk and sample type
        new_labels = format_xtick_labels(ax.get_xticklabels())
        ax.set_xticklabels(new_labels)

        # add the tick locator:
        if molecule == "aa":
            if vrc01 == "vrc01":
                ax.set_ylim(-0.5, 16)
                ax.yaxis.set_major_locator(plt.MultipleLocator(2))
                ax.yaxis.set_minor_locator(plt.MultipleLocator(1))
                ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: str(int(y))))
            else:
                ax.set_ylim(-0.5, 28)
                ax.yaxis.set_major_locator(plt.MultipleLocator(4))
                ax.yaxis.set_minor_locator(plt.MultipleLocator(2))
            ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: str(int(y)) if y != 28 else "    28"))
        else:
            # if vrc01 == "vrc01":
            # ax.set_ylim(0, 6.5)
            # else:
            # ax.set_ylim(0, 23)
            ax.yaxis.set_major_locator(plt.MultipleLocator(2))
            ax.yaxis.set_minor_locator(plt.MultipleLocator(1))

        # set column tittles that are dosage groups
        if ax_index == 0:
            ax.set_title(r"20 $\mu$g Dose", pad=20, fontsize=12)
        if ax_index == 1:
            ax.set_title(r"100 $\mu$g Dose", pad=20, fontsize=12)

    # finally remove spnes, not sure if this is necesary, but it is necessary to remove all legends
    for _, ax in enumerate(concatenate((row_ax1, row_ax2))):
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.get_legend().remove()

    # add a supfigure title slightly off centered
    fig.suptitle(sup_title, ha="center", x=0.53)
    fig.tight_layout()
    return fig
