# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
from __future__ import annotations
from itertools import product

# stdlib
import logging
from pathlib import Path

# third party
import rich_click as click
import seaborn as sns
import numpy as np
from g001.figures.binding import plot_boosting_binding, plot_gt8_binding
from g001.figures.features import plot_sequence_features
from g001.figures.polyclonality import plot_class_clonality

# G001 package
from g001.data import Data
from g001.figures.frequency import (
    plot_count_frequency_response,
    plot_flow_frequencies,
    plot_vrc01_among_epitope_specific,
)
from g001.figures.overview import plot_overview
from g001.figures.somatic import plot_somatic_mutation_frequencies, plot_somatic_mutation_frequencies_violin
from g001.sequence_pipeline.main import run_sequence_analysis
from g001.utils.latex_modifiers import make_tables
from g001.utils.rscript_wrapper import RScript


@click.group()
@click.pass_context
@click.option(
    "--data-path", "-d", type=click.Path(exists=True, dir_okay=True), default=Path("./data"), show_default=True
)
def main(ctx: click.Context, data_path: Path) -> None:
    """
    Generate G001
    """
    logging.basicConfig(
        level=logging.INFO,
        datefmt="%m-%d-%Y %I:%M%p",
        format="%(name)s:%(levelname)s:%(asctime)s--> %(message)s",
    )
    data = Data(data_path)

    # make empty context object
    ctx.obj = {"data": data}
    pass


@main.command("sequence_analysis")
@click.pass_context
@click.option(
    "--data-path",
    "-d",
    type=click.Path(dir_okay=True, readable=True, exists=True, resolve_path=True),
    required=True,
    help="Path to the sequence data directory",
)
@click.option(
    "--outpath",
    "-o",
    type=click.Path(dir_okay=True, readable=True, resolve_path=True),
    required=True,
    help="Outpath to sequence_analysis",
    default="sequence_analysis_output",
)
@click.option("--resume", "-r", is_flag=True, help="""Resume pipeline from .cache path""")
@click.option("--skip-csv", is_flag=True, help="""Skip copying feathers to .csv in final path""")
def click_run_sa(ctx: click.Context, data_path: str | Path, outpath: Path, resume: bool, skip_csv: bool) -> None:
    """
    Run sequence analysis
    """
    logger = logging.getLogger("G001")
    data = Data(Path(data_path))
    logger.info(
        f"Running complete sequence analysis pipeline using data: {data.sequence_data_paths.base_sequence_path}"
    )

    # create output directory and cache path
    outpath = Path(outpath)
    cache_path = outpath / Path(".cache")

    # if they dont exist, create them
    if not outpath.exists():
        logger.info(f"Creating output directory: {outpath}")
        outpath.mkdir(parents=True)
    if not cache_path.exists():
        logger.info(f"Creating cache directory: {cache_path}")
        cache_path.mkdir(parents=True)

    # run sequene analysis pipeline
    run_sequence_analysis(data, cache_path, resume, skip_csv)


@main.command("process-flow")
@click.pass_context
@click.option(
    "--verbose",
    "-v",
    is_flag=True,
    required=False,
    default=False,
    help="Prints logs to console",
)
@click.option(
    "--site",
    "-s",
    required=True,
    type=click.Choice(["FHCRC", "VRC"]),
)
@click.option(
    "--manifest",
    "-m",
    type=click.Path(dir_okay=False, readable=True, exists=True, resolve_path=True),
    required=True,
    show_default=True,
    help="Path to read the flow manifest for specified site",
)
@click.option(
    "--flow-input-dir",
    "-i",
    type=click.Path(dir_okay=True, readable=True, exists=True, resolve_path=True),
    required=False,
    default="flow_input",
    show_default=True,
    help="Path to read the flow_input directory",
)
@click.option(
    "--flow-output-dir",
    "-o",
    type=click.Path(dir_okay=True, readable=True, resolve_path=True),
    required=False,
    default=Path("flow_output"),
    show_default=True,
    help="Path to create the flow_output directory",
)
@click.option("-f", "--force_overwrite_output_dir", default=True, is_flag=True, help="Force overwrite output directory")
def process_flow(
    ctx: click.Context,
    verbose: bool,
    site: str,
    manifest: Path | str,
    flow_input_dir: Path | str,
    flow_output_dir: Path | str,
    force_overwrite_output_dir: str,
) -> None:
    """
    Process flow data
    """
    if not Path(flow_output_dir).exists():
        Path(flow_output_dir).mkdir()
    RScript(verbose=verbose).flow_processing(
        site=site,
        manifest=manifest,
        flow_input_dir=flow_input_dir,
        flow_output_dir=flow_output_dir,
        force_overwrite_output_dir=force_overwrite_output_dir,
    )


@main.command("collate")
@click.pass_context
@click.option(
    "--verbose",
    "-v",
    is_flag=True,
    required=False,
    default=False,
    help="Prints final command that was run to console",
)
@click.option(
    "--fhcrc-manifest",
    "-1",
    type=click.Path(dir_okay=False, readable=True, exists=True, resolve_path=True),
    required=False,
    show_default=True,
    help="Path to read the flow manifest for specified site",
)
@click.option(
    "--vrc-manifest",
    "-2",
    type=click.Path(dir_okay=False, readable=True, exists=True, resolve_path=True),
    required=False,
    show_default=True,
    help="Path to read the flow manifest for specified site",
)
@click.option(
    "--flow-processed-dir",
    "-f",
    type=click.Path(dir_okay=True, readable=True, exists=True, resolve_path=False),
    required=False,
    default=None,
    help="Path to the processed flow directory containing both FHCRC and VRC data",
)
@click.option(
    "--collated-output-dir",
    "-o",
    type=click.Path(dir_okay=True, readable=True, exists=False, resolve_path=False),
    required=False,
    default="collated_flow",
    help="Output path for the collated flow data",
)
def collate(
    ctx: click.Context,
    verbose: bool,
    fhcrc_manifest: Path | str,
    vrc_manifest: Path | str,
    flow_processed_dir: Path | str,
    collated_output_dir: Path | str,
) -> None:
    """Collated Flow Data from sites FHCRC & VRC"""
    data = ctx.obj["data"]
    if not flow_processed_dir:
        flow_processed_dir = data.get_processed_flow_paths()
    if not fhcrc_manifest:
        fhcrc_manifest = flow_processed_dir / Path("fhrc/fhcrc_manifest.csv")
    if not vrc_manifest:
        vrc_manifest = flow_processed_dir / Path("vrc/vrc_manifest.csv")
    if not Path(collated_output_dir).exists():
        Path(collated_output_dir).mkdir()
    RScript(verbose=verbose).collate_flow(
        fhcrc_manifest=Path(fhcrc_manifest),
        vrc_manifest=Path(vrc_manifest),
        flow_processed_dir=Path(flow_processed_dir),
        collated_output_dir=Path(collated_output_dir),
    )


@main.command("combine")
@click.pass_context
@click.option(
    "--verbose",
    "-v",
    is_flag=True,
    required=False,
    default=False,
    help="Prints final command that was run to console",
)
@click.option(
    "--fhcrc-manifest",
    "-1",
    type=click.Path(dir_okay=False, readable=True, exists=True, resolve_path=True),
    required=False,
    show_default=True,
    help="Path to read the flow manifest for specified site",
)
@click.option(
    "--vrc-manifest",
    "-2",
    type=click.Path(dir_okay=False, readable=True, exists=True, resolve_path=True),
    required=False,
    show_default=True,
    help="Path to read the flow manifest for specified site",
)
@click.option(
    "--seq-dir",
    "-s",
    type=click.Path(dir_okay=True, readable=True, exists=True, resolve_path=False),
    required=False,
    default=None,
    help="Path to the processed flow directory containing both FHCRC and VRC data",
)
@click.option(
    "--collated-output-dir",
    "-c",
    type=click.Path(dir_okay=True, readable=True, exists=False, resolve_path=False),
    required=False,
    default="collated_flow",
    help="Output path for the collated flow data",
)
@click.option(
    "--combined-output-dir",
    "-o",
    type=click.Path(dir_okay=True, readable=True, exists=False, resolve_path=False),
    required=False,
    default="combined_flow_seq",
    help="Output path containing the combined flow data and sequence data",
)
def collate(
    ctx: click.Context,
    verbose: bool,
    fhcrc_manifest: Path,
    vrc_manifest: Path,
    seq_dir: Path,
    collated_output_dir: Path,
    combined_output_dir: Path,
) -> None:
    """Collated Flow Data from sites FHCRC & VRC"""
    data = ctx.obj["data"]
    flow_processed_dir = data.get_processed_flow_paths()
    if not fhcrc_manifest:
        fhcrc_manifest = flow_processed_dir / Path("fhrc/fhcrc_manifest.csv")
    if not vrc_manifest:
        vrc_manifest = flow_processed_dir / Path("vrc/vrc_manifest.csv")
    if not seq_dir:
        seq_dir = data.sequence_data_paths.base_sequence_path
    if not Path(combined_output_dir).exists():
        Path(combined_output_dir).mkdir()
    RScript(verbose=verbose).combine_data_seq(
        fhcrc_manifest=Path(fhcrc_manifest),
        vrc_manifest=Path(vrc_manifest),
        seq_dir=Path(seq_dir),
        collated_output_dir=Path(collated_output_dir),
        combined_output_dir=Path(combined_output_dir),
    )


@main.command("combine")
@click.pass_context
@click.option(
    "--verbose",
    "-v",
    is_flag=True,
    required=False,
    default=False,
    help="Prints final command that was run to console",
)
@click.option(
    "--fhcrc-manifest",
    "-1",
    type=click.Path(dir_okay=False, readable=True, exists=True, resolve_path=True),
    required=False,
    show_default=True,
    help="Path to read the flow manifest for specified site",
)
@click.option(
    "--vrc-manifest",
    "-2",
    type=click.Path(dir_okay=False, readable=True, exists=True, resolve_path=True),
    required=False,
    show_default=True,
    help="Path to read the flow manifest for specified site",
)
@click.option(
    "--sequence-dir",
    "-s",
    type=click.Path(dir_okay=True, readable=True, exists=True, resolve_path=True),
    required=False,
    show_default=True,
    help="Path to read the sequencing data",
)
@click.option(
    "--collated-dir",
    "-c",
    type=click.Path(dir_okay=True, readable=True, exists=True, resolve_path=True),
    required=False,
    show_default=True,
    help="path to the collated directory",
)
@click.option(
    "--combined-output-dir",
    "-o",
    type=click.Path(dir_okay=True, readable=True, exists=False, resolve_path=True),
    required=False,
    default="combined_data",
)
def combine(
    ctx: click.Context,
    verbose: bool,
    fhcrc_manifest: Path,
    vrc_manifest: Path,
    sequence_dir: Path,
    collated_dir: Path,
    combined_output_dir: Path,
) -> None:
    """Collated Flow Data from sites FHCRC & VRC"""
    data: Data = ctx.obj["data"]
    flow_processed_dir = data.get_processed_flow_paths()
    if not fhcrc_manifest:
        fhcrc_manifest = flow_processed_dir / Path("fhrc/fhcrc_manifest.csv")
    if not vrc_manifest:
        vrc_manifest = flow_processed_dir / Path("vrc/vrc_manifest.csv")
    if not sequence_dir:
        sequence_dir = data.get_data_sequence_path()
    if not collated_dir:
        collated_dir = data.get_collated_data_path()
    if not Path(combined_output_dir).exists():
        Path(combined_output_dir).mkdir()
    RScript(verbose=verbose).combine_flow_and_sequence(
        fhcrc_manifest=Path(fhcrc_manifest),
        vrc_manifest=Path(vrc_manifest),
        sequence_dir=Path(sequence_dir),
        collated_dir=Path(collated_dir),
        combined_output_dir=Path(combined_output_dir),
    )


@main.group()
def figures():
    """
    Generate for figures.
    """
    sns.set_context("paper", font_scale=1)
    sns.set_style("ticks")
    sns.set_style({"font.family": "Arial"})
    import warnings

    warnings.simplefilter("ignore")
    np.random.seed(1000)
    pass


@figures.command("fig1")
@click.pass_context
@click.option(
    "--outpath",
    "-o",
    default="figures/figure1",
    type=click.Path(file_okay=True, resolve_path=True),
    help="Output path for figure",
    show_default=True,
)
def plot_figure_1(ctx: click.Context, outpath: str | Path) -> None:
    data = ctx.obj["data"]
    figure = plot_flow_frequencies(data)
    figure.savefig(outpath + ".png", dpi=300)
    click.echo(f"Figure 1 generated to {outpath}.png")


@figures.command("fig2")
@click.pass_context
@click.option(
    "--outpath",
    "-o",
    default="figures/figure2",
    type=click.Path(file_okay=True, resolve_path=True),
    help="Output path for figure",
    show_default=True,
)
def plot_figure_2(ctx: click.Context, outpath: str | Path) -> None:
    data = ctx.obj["data"]
    figure = plot_count_frequency_response(data)
    figure.savefig(outpath + ".png", dpi=300)
    click.echo(f"Figure 2A-C generated to {outpath}.png")

    figure_d = plot_overview(data)
    figure_d.savefig(outpath + "_d.png", dpi=300)
    click.echo(f"Figure 2D generated to {outpath}_d.png")


@figures.command("fig3")
@click.pass_context
@click.option(
    "--outpath",
    "-o",
    default="figures/figure3",
    type=click.Path(file_okay=True, resolve_path=True),
    help="Output path for figure",
    show_default=True,
)
def plot_figure_3(ctx: click.Context, outpath: str | Path) -> None:
    data = ctx.obj["data"]
    figure = plot_class_clonality(data, True)
    figure.savefig(outpath + ".png", dpi=300)
    click.echo(f"Figure 3 generated to {outpath}.png")


@figures.command("fig4")
@click.pass_context
@click.option(
    "--outpath",
    "-o",
    default="figures/figure4",
    type=click.Path(file_okay=True, resolve_path=True),
    help="Output path for figure",
    show_default=True,
)
def plot_figure_4(ctx: click.Context, outpath: str | Path) -> None:
    data = ctx.obj["data"]
    figure = plot_vrc01_among_epitope_specific(data)
    figure.savefig(outpath + ".png", dpi=300)
    click.echo(f"Figure 4 generated to {outpath}.png")


@figures.command("fig5")
@click.pass_context
@click.option(
    "--aa-outpath",
    "-o",
    default="figures/figure5",
    type=click.Path(file_okay=True, resolve_path=True),
    help="Output path for figure",
    show_default=True,
)
@click.option(
    "--nt-outpath",
    "-n",
    default="figures/figureS20",
    type=click.Path(file_okay=True, resolve_path=True),
    help="Output path for figure S5",
    show_default=True,
)
def plot_figure_5_and_s26(ctx: click.Context, aa_outpath: str, nt_outpath: str) -> None:
    data: Data = ctx.obj["data"]
    plot_params = data.plot_parameters
    unblind_df = data.get_unblinded_sequences()
    pallete = plot_params.get_pallete()
    vrc01_dose_palette = {
        "Low Dose": pallete["Low Dose"]["vrc01_class"],
        "High Dose": pallete["High Dose"]["vrc01_class"],
    }
    nonvrc01_dose_palette = {
        "Low Dose": pallete["Low Dose"]["nonvrc01_class"],
        "High Dose": pallete["High Dose"]["nonvrc01_class"],
    }

    for ab_type, molecule in list(product(["vrc01", "nonvrc01"], ["nt", "aa"])):
        if ab_type == "vrc01":
            sup_title = "VRC01-class BCRs"
            palette = vrc01_dose_palette
        else:
            sup_title = "non-VRC01-class BCRs"
            palette = nonvrc01_dose_palette
        if molecule == "nt":
            outpath_ = nt_outpath
        else:
            outpath_ = aa_outpath
        fig = plot_somatic_mutation_frequencies(unblind_df, palette, ab_type, molecule, sup_title)
        new_outpath = outpath_ + f"_{ab_type}_{molecule}.png"
        fig.savefig(new_outpath, dpi=300)
        if molecule == "nt":
            _figure = "S20"
        else:
            _figure = "2"
        click.echo(f"Figure {_figure} generated to {new_outpath}")
        fig = plot_somatic_mutation_frequencies_violin(unblind_df, palette, ab_type, molecule, sup_title)
        new_outpath = outpath_ + f"_{ab_type}_{molecule}_violin.png"
        click.echo(f"Figure {_figure} generated to {new_outpath}")
        fig.savefig(new_outpath, dpi=300)


@figures.command("fig6")
@click.pass_context
@click.option(
    "--outpath",
    "-o",
    default="figures/figure6",
    type=click.Path(file_okay=True, resolve_path=True),
    help="Output path for figure",
    show_default=True,
)
def plot_figure_6(ctx: click.Context, outpath: str | Path) -> None:
    data = ctx.obj["data"]
    figure = plot_sequence_features(data)
    figure.savefig(outpath + ".png", dpi=300)
    click.echo(f"Figure 6 generated to {outpath}.png")


@figures.command("fig7")
@click.pass_context
@click.option(
    "--outpath",
    "-o",
    default="figures/figure7",
    type=click.Path(file_okay=True, resolve_path=True),
    help="Output path for figure",
    show_default=True,
)
def plot_figure_6(ctx: click.Context, outpath: str | Path) -> None:
    data = ctx.obj["data"]
    figure = plot_gt8_binding(data)
    figure.savefig(outpath + ".png", dpi=300)
    click.echo(f"Figure 7 generated to {outpath}.png")


@figures.command("fig8")
@click.pass_context
@click.option(
    "--outpath",
    "-o",
    default="figures/figure8",
    type=click.Path(file_okay=True, resolve_path=True),
    help="Output path for figure",
    show_default=True,
)
def plot_figure_8(ctx: click.Context, outpath: str | Path) -> None:
    data = ctx.obj["data"]
    figure = plot_boosting_binding(data)
    figure.savefig(outpath + ".png", dpi=300)
    click.echo(f"Figure 8 generated to {outpath}.png")


@main.command("supptables")
@click.pass_context
@click.option(
    "--outpath",
    "-o",
    type=click.Path(file_okay=True, resolve_path=True),
    help="Output path for Sup Tables for Paper 1",
    show_default=True,
    required=True,
)
@click.option("-c", "--cleanheaders", is_flag=True, default=False)
def generate_st1(ctx: click.Context, cleanheaders: bool, outpath: str) -> None:
    data = ctx.obj["data"]
    if not outpath:
        table_output = data.table_data_paths.table_out
    else:
        if not Path(outpath).exists():
            Path(outpath).mkdir(parents=True)
            (Path(outpath) / Path("pdfs")).mkdir()
        table_output = Path(outpath)
    table_src = data.table_data_paths.table_src

    pdf_order = [
        # Static files with Word file src
        ("S1", "", "misc/demographics.docx"),
        ("S2", "", "misc/G001 Schedule of Procedures from Protocol from Vince 01Aug2021.docx"),
        ("S3", "", "misc/G001 Trials.pdf"),  # made from df
        ("S4", "", "misc/solicited symptoms.docx"),
        ("S5", "", "misc/unsolicited AEs.docx"),
        ("S6", "", "misc/unsolicited AEs related to study procedures.docx"),
        ("S7", "", "misc/response-rates-IQR-5-antigens-epitopes.pdf"),
        ("S8", "", "misc/AUTC-mag-weeks2_10-5-antigens-epitopes.pdf"),
        ("S9", "", "misc/Response rate summary.tex"),
        # Kristen tables
        ("S10", "", "kristen-tables/Memory B Cell and Plasmablast Flow Cytometry Panel.pdf"),
        ("S11", "", "kristen-tables/Germinal Center B Cell LN FNA Flow Cytometry Panel.pdf"),
        ("S12", "", "kristen-tables/Table of Reagents.pdf"),
        ("S13", "", "kristen-tables/PCR-IgH Primary (Round 1).pdf"),
        ("S14", "", "kristen-tables/PCR-IgH Nested (Round 2).pdf"),
        ("S15", "", "kristen-tables/PCR-Igκ Primary (Round 1).pdf"),
        ("S16", "", "kristen-tables/PCR-Igκ Nested (Round 2).pdf"),
        ("S17", "", "kristen-tables/PCR-Igλ Primary (Round 1).pdf"),
        ("S18", "", "kristen-tables/PCR-Igλ Nested (Round 2).pdf"),
        ("S19", "", "kristen-tables/PCR constructs.pdf"),
        # Files built from Cell Report
        ("S20", "B cell Report Table 2", "individual_tables/bcell-tab-02.tex"),
        ("S21", "B cell Report Table 30", "individual_tables/bcell-tab-30.tex"),
        ("S22", "B cell Report Table 18", "individual_tables/bcell-tab-18.tex"),
        ("S23", "B cell Report Table 29", "individual_tables/bcell-tab-29.tex"),
        ("S24", "B cell Report Table 17", "individual_tables/bcell-tab-17.tex"),
        ("S25", "B cell Report Table 31", "individual_tables/bcell-tab-31.tex"),
        ("S26", "B cell Report Table 32", "individual_tables/bcell-tab-32.tex"),
        ("S27", "B cell Report Table 19", "individual_tables/bcell-tab-19.tex"),
        ("S28", "B cell Report Table 20", "individual_tables/bcell-tab-20.tex"),
        ("S29", "B cell Report Table 33", "individual_tables/bcell-tab-33.tex"),
        ("S30", "B cell Report Table 21", "individual_tables/bcell-tab-21.tex"),
        ("S31", "B cell Report Table 44", "individual_tables/bcell-tab-44.tex"),
        ("S32", "B cell Report Table 5", "individual_tables/bcell-tab-57.tex"),
        ("S33", "B cell Report Table 48", "individual_tables/bcell-tab-58.tex"),
        ("S34", "B cell Report Table 50", "individual_tables/bcell-tab-59.tex"),
        ("S35", "B cell Report Table 6", "individual_tables/bcell-tab-06.tex"),
        ("S36", "B cell Report Table 5", "individual_tables/bcell-tab-05.tex"),
        ("S37", "B cell Report Table 48", "individual_tables/bcell-tab-48.tex"),
        ("S38", "B cell Report Table 50", "individual_tables/bcell-tab-50.tex"),
        ("S39", "B cell Report Table 47", "individual_tables/bcell-tab-47.tex"),
        ("S40", "B cell Report Table 52", "individual_tables/bcell-tab-52.tex"),
        ("S41", "B cell Report Table 49", "individual_tables/bcell-tab-49.tex"),
        ("S42", "B cell Report Table 12", "individual_tables/bcell-tab-12.tex"),
        # Files built from Genotype Report
        ("S43", "Genotype Report Table", "individual_tables/genotype-tab-01_02.tex"),
        # Merged from Paper 2
        ("S44", "Percent Mutation Report Table 2", "individual_tables/percent-mut-tab-02.tex"),
        ("S45", "Percent Mutation Report Table 8", "individual_tables/percent-mut-tab-08.tex"),
        ("S46", "Percent Mutation Report Table 4", "individual_tables/percent-mut-tab-04.tex"),
        ("S47", "Percent Mutation Report Table 5", "individual_tables/percent-mut-tab-05.tex"),
        ("S48", "Percent Mutation Report Table 6", "individual_tables/percent-mut-tab-06.tex"),
        ("S49", "Percent Mutation Report Table 7", "individual_tables/percent-mut-tab-07.tex"),
        ("S50", "", "misc/BCRs detected.pdf"),  # static
    ]
    make_tables(table_src, table_output, pdf_order, cleanheaders=cleanheaders, paper_numnber=1)


if __name__ == "__main__":
    main()
