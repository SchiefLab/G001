# -*- coding: utf-8 -*-
import logomaker
import pandas as pd
from matplotlib import ticker as mtick


def plot_sequence_logo(
    ax, show_left_spine=False, ytitle=None, title=None, char_set=["ABCDEF"], xticks=[], ytitle_font=8
):
    color_scheme = "NajafabadiEtAl2017"

    # make matrix for character set
    mat_df = logomaker.alignment_to_matrix(char_set)
    mat_df = logomaker.transform_matrix(mat_df, normalize_values=True)
    # make and style logo
    logomaker.Logo(mat_df, ax=ax, color_scheme=color_scheme, show_spines=False)
    if xticks:
        ax.xaxis.set_major_locator(mtick.FixedLocator(range(0, 5)))
    else:
        ax.set_xticks([])
    ax.set_yticks([])
    ax.set_yticklabels([])
    ax.set_xticklabels(xticks, rotation=0, fontsize=7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if show_left_spine:
        ax.spines["left"].set_visible(True)
    if title:
        ax.set_title(title, fontsize=10)
    if ytitle:
        ax.set_ylabel(ytitle, fontsize=ytitle_font, rotation=0, ha="right", va="center")


def get_sequences_from_df(df, is_vrc01_class, dose_group, locus):
    _seqs = (
        df.query("vaccine_group=='vaccine'")
        .query(f"dose_group=='{dose_group}'")
        .query(f"is_vrc01_class=={is_vrc01_class}")
        .query("has_5_len_lcdr3")
        .query(f"locus_light=='{locus}'")["cdr3_aa_light"]
    ).to_list()
    return _seqs


def plot_seq_logo_panel(
    fig,
    gs,
    kappa_vrc01_class_df: pd.DataFrame,
    lambda_vrc01_class_df: pd.DataFrame,
    oas_5_len_df: pd.DataFrame,
    plottable_df: pd.DataFrame,
    human_naive_5_len_df: pd.DataFrame,
):
    """seq logo madness"""
    seq_logo_panel_gs = gs

    low_dose_label = r"20 $\mu$g"
    high_dose_label = r"100 $\mu$g"

    # establish axis for vrc01 class
    vrc01_class_title_axis = fig.add_subplot(seq_logo_panel_gs[0, :])
    bnabs_kappa_vrc01_class_axis = fig.add_subplot(seq_logo_panel_gs[1, 0])
    bnabs_lambda_vrc01_class_axis = fig.add_subplot(seq_logo_panel_gs[1, 1])
    low_dose_kappa_vrc01_class_axis = fig.add_subplot(seq_logo_panel_gs[2, 0])
    low_dose_lambda_vrc01_class_axis = fig.add_subplot(seq_logo_panel_gs[2, 1])
    high_dose_kappa_vrc01_class_axis = fig.add_subplot(seq_logo_panel_gs[3, 0])
    high_dose_lambda_vrc01_class_axis = fig.add_subplot(seq_logo_panel_gs[3, 1])
    humain_naive_kappa_axis = fig.add_subplot(seq_logo_panel_gs[4, 0])
    humain_naive_lambda_axis = fig.add_subplot(seq_logo_panel_gs[4, 1])

    # nonvrc01 class seq logos
    buffer = fig.add_subplot(seq_logo_panel_gs[5, :])
    non_vrc01_class_title_axis = fig.add_subplot(seq_logo_panel_gs[6, :])
    bnabs_kappa_vrc01_class_axis_2 = fig.add_subplot(seq_logo_panel_gs[7, 0])
    bnabs_lambda_vrc01_class_axis_2 = fig.add_subplot(seq_logo_panel_gs[7, 1])
    low_dose_kappa_nonvrc01_class_axis = fig.add_subplot(seq_logo_panel_gs[8, 0])
    low_dose_lambda_nonvrc01_class_axis = fig.add_subplot(seq_logo_panel_gs[8, 1])
    high_dose_kappa_nonvrc01_class_axis = fig.add_subplot(seq_logo_panel_gs[9, 0])
    high_dose_lambda_nonvrc01_class_axis = fig.add_subplot(seq_logo_panel_gs[9, 1])
    oas_5_len_kappa_axis = fig.add_subplot(seq_logo_panel_gs[10, 0])
    oas_5_len_lambda_axis = fig.add_subplot(seq_logo_panel_gs[10, 1])

    # Collect all sequences we need
    kappa_vrc01_class_seqs = kappa_vrc01_class_df["cdr3_aa_no_gaps"].to_list()
    lambda_vrc01_class_seqs = lambda_vrc01_class_df["cdr3_aa_no_gaps"].to_list()

    # human naive - this is linked dataframe so need cdr3_aa_light
    human_naive_5_len_kappa = human_naive_5_len_df.query("locus_light=='IGK'")["cdr3_aa_light"].to_list()
    human_naive_5_len_lambda = human_naive_5_len_df.query("locus_light=='IGL'")["cdr3_aa_light"].to_list()

    # oas 5 len
    oas_5_len_kappa_seqs = oas_5_len_df.query("locus=='IGK'")["cdr3_aa"].to_list()
    oas_5_len_lambda_seqs = oas_5_len_df.query("locus=='IGL'")["cdr3_aa"].to_list()

    # low dose
    low_dose_kappa_vrc01 = get_sequences_from_df(plottable_df, True, "Low Dose", "IGK")
    low_dose_lambda_vrc01 = get_sequences_from_df(plottable_df, True, "Low Dose", "IGL")

    # low dose non vrc01
    low_dose_kappa_nonvrc01 = get_sequences_from_df(plottable_df, False, "Low Dose", "IGK")
    low_dose_lambda_nonvrc01 = get_sequences_from_df(plottable_df, False, "Low Dose", "IGL")

    # high dose
    high_dose_kappa_vrc01 = get_sequences_from_df(plottable_df, True, "High Dose", "IGK")
    high_dose_lambda_vrc01 = get_sequences_from_df(plottable_df, True, "High Dose", "IGL")

    # high dose non vrc01
    high_dose_kappa_nonvrc01 = get_sequences_from_df(plottable_df, False, "High Dose", "IGK")
    high_dose_lambda_nonvrc01 = get_sequences_from_df(plottable_df, False, "High Dose", "IGL")

    vrc01_class_title_axis.annotate(
        "VRC01 Class",
        xy=(0.5, 1),
        xytext=(0.5, 1.3),
        xycoords="axes fraction",
        ha="center",
        va="top",
        fontsize=10,
    )
    vrc01_class_title_axis.axis("off")
    non_vrc01_class_title_axis.annotate(
        "Non-VRC01-class 5aa LCDR3s",
        xy=(0.5, 1.0),
        xytext=(0.5, 1.3),
        xycoords="axes fraction",
        ha="center",
        va="top",
        fontsize=10,
    )
    non_vrc01_class_title_axis.axis("off")
    buffer.axis("off")

    # bnabs kappa
    plot_sequence_logo(
        bnabs_kappa_vrc01_class_axis,
        show_left_spine=True,
        char_set=kappa_vrc01_class_seqs,
        title="Kappa",
        ytitle="bnAbs",
    )

    # bnabs lambda
    plot_sequence_logo(
        bnabs_lambda_vrc01_class_axis,
        show_left_spine=False,
        char_set=lambda_vrc01_class_seqs,
        title="Lambda",
    )

    # vrc01 low dose kappa
    plot_sequence_logo(
        low_dose_kappa_vrc01_class_axis,
        show_left_spine=True,
        char_set=low_dose_kappa_vrc01,
        ytitle=low_dose_label,
    )

    # vrc01 low dose lambda
    plot_sequence_logo(
        low_dose_lambda_vrc01_class_axis,
        show_left_spine=False,
        char_set=low_dose_lambda_vrc01,
    )

    # vrc01 high dose kappa
    plot_sequence_logo(
        high_dose_kappa_vrc01_class_axis,
        show_left_spine=True,
        char_set=high_dose_kappa_vrc01,
        ytitle=high_dose_label,
    )

    # vrc01 high dose lambdas
    plot_sequence_logo(
        high_dose_lambda_vrc01_class_axis,
        show_left_spine=False,
        char_set=high_dose_lambda_vrc01,
    )

    # plot vrc01 control
    plot_sequence_logo(
        humain_naive_kappa_axis,
        xticks=list(range(93, 98)),
        show_left_spine=True,
        char_set=human_naive_5_len_kappa,
        ytitle="Human naive\nGT8 binders",
    )

    humain_naive_kappa_axis.annotate(
        "Light Chain Position",
        xy=(1, -0.5),
        xytext=(1.05, -0.4),
        xycoords="axes fraction",
        ha="center",
        va="top",
        fontsize=10,
    )

    plot_sequence_logo(
        humain_naive_lambda_axis,
        xticks=list(range(93, 98)),
        show_left_spine=False,
        char_set=human_naive_5_len_lambda,
    )

    # do non vrc01
    plot_sequence_logo(
        bnabs_kappa_vrc01_class_axis_2,
        show_left_spine=True,
        char_set=kappa_vrc01_class_seqs,
        title="Kappa",
        ytitle="VRC01-class\nbnAbs",
    )

    # bnabs lambda
    plot_sequence_logo(
        bnabs_lambda_vrc01_class_axis_2,
        show_left_spine=False,
        char_set=lambda_vrc01_class_seqs,
        title="Lambda",
    )
    plot_sequence_logo(
        low_dose_kappa_nonvrc01_class_axis,
        show_left_spine=True,
        char_set=low_dose_kappa_nonvrc01,
        ytitle=low_dose_label,
    )
    plot_sequence_logo(
        low_dose_lambda_nonvrc01_class_axis,
        show_left_spine=False,
        char_set=low_dose_lambda_nonvrc01,
    )
    plot_sequence_logo(
        high_dose_kappa_nonvrc01_class_axis,
        show_left_spine=True,
        char_set=high_dose_kappa_nonvrc01,
        ytitle=high_dose_label,
    )
    plot_sequence_logo(
        high_dose_lambda_nonvrc01_class_axis,
        show_left_spine=False,
        char_set=high_dose_lambda_nonvrc01,
    )
    plot_sequence_logo(
        oas_5_len_kappa_axis,
        show_left_spine=True,
        char_set=oas_5_len_kappa_seqs,
        xticks=list(range(93, 98)),
        ytitle="OAS 5aa L3",
    )
    plot_sequence_logo(
        oas_5_len_lambda_axis,
        show_left_spine=False,
        char_set=oas_5_len_lambda_seqs,
        xticks=list(range(93, 98)),
    )
    oas_5_len_kappa_axis.annotate(
        "Light Chain Position",
        xy=(1, -0.5),
        xytext=(1.05, -0.4),
        xycoords="axes fraction",
        ha="center",
        va="top",
        fontsize=10,
    )
