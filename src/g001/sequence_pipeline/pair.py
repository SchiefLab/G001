from distutils.spawn import find_executable
from g001.data import Data
import pandas as pd
import logging
from typing import Tuple

logger = logging.getLogger("Pair")


def pair_chains(unpaired_dataframe: pd.DataFrame) -> pd.DataFrame:
    """Given the unpaired dataframe, pair heavy and light chains

    Parameters
    ----------
    unpaired_dataframe : pd.DataFrame
        The unpaired dataframe

    Returns
    -------
    pd.DataFrame
        The paired dataframe
    """
    # dropped_for na are the sequences that are IGH_IGK or IGH_IGL
    pc = unpaired_dataframe[unpaired_dataframe["dropped_for"].isna()].copy()

    logger.info(f"{len(pc)} sequences are IGH_IGK or IGH_IGL...paring")
    heavy = pc.query("locus=='IGH'")
    light = pc.query("locus=='IGK' or locus=='IGL'")
    group = ["pub_id", "timepoint", "plate", "well", "plate_type", "tissue"]

    # dropping the dropped for column since they are all na
    clean_pairs = heavy.drop("dropped_for", axis=1).merge(
        light.drop("dropped_for", axis=1), on=group, suffixes=["_heavy", "_light"]
    )
    logger.info(f"Paired: {len(clean_pairs):,}")
    return clean_pairs


def find_pairing_type(remaining_dataframe: pd.DataFrame) -> pd.DataFrame:
    """Look into each well and find the pairing type. It could be multiple heavy...mutli light, heavy_kappa...etc.
    Will add a field called pairing type

    Parameters
    ----------
    remaining_dataframe : pd.DataFrame
        The pairing candidate dataframe
    Returns
    -------
    pd.DataFrame
        The same dataframe with paring type added
    """

    # turn locus into a category so we can sort it
    remaining_dataframe["locus"] = pd.Categorical(
        remaining_dataframe["locus"], categories=["IGH", "IGK", "IGL"], ordered=True
    )

    # for the index, reset it so its a column...will need later
    remaining_dataframe = remaining_dataframe.reset_index()

    # for each group at a particular well, we will join the locus by a "_". So you will get IGH_IGK, IGH_IGL...etc
    pairing_candidates = (
        remaining_dataframe.sort_values("locus")
        .groupby(["pub_id", "timepoint", "plate", "well"])
        .apply(lambda x: pd.Series({"pairing_type": "_".join(x["locus"].to_list())}))
    )

    # reset the group and merge back on remaining dataframe...
    pairing_candidates_merged = remaining_dataframe.merge(
        pairing_candidates.reset_index(), on=["pub_id", "timepoint", "plate", "well"]
    )

    pairing_type_series = pairing_candidates_merged["pairing_type"].value_counts()
    logger.info(f"Pairing type results:\n{pairing_type_series}")

    # because we need to fill in the unpaired dataframe, we need to set the column index back to the actual index
    pairing_candidates_merged = pairing_candidates_merged.set_index("index")
    pairing_candidates_merged.index.name = ""
    return pairing_candidates_merged


def drop_bad_samples(working_dataframe: pd.DataFrame) -> pd.DataFrame:
    """This function does two things...

    1. It will add a field called "dropped for". This field tells the final fate of each sequence and why it was dropped.
    2. For any sequences that are not dropped, it will determine the pairing type. IGH_IGK, IGH_IGL. For anything that is not that pairing type, it will update the sample to tell it that the sample is "unpaired"

    Parameters
    ----------
    working_dataframe : pd.DataFrame
        the working dataframe with all the sequences

    Returns
    -------
    pd.DataFrame
        The working dataframe with a "dropped for" field added with the final fates of all the sequences
    """
    remaining_dataframe = working_dataframe.copy()
    unpaired_dataframe = working_dataframe.copy()

    # Step 0 - Only accept latest reps per our discussion
    grouping_key = ["pub_id", "timepoint", "plate", "well", "chain"]

    # get the indexes of the last replicate
    last_rep_index = unpaired_dataframe.sort_values("replicate").groupby(grouping_key).tail(1).index  # type: ignore

    # get the indexes that are not last rep
    not_last_rep = unpaired_dataframe.index.difference(last_rep_index)
    # logger.info(f"Saving {len(not_last_rep.intersection(reseq_corrected_index))}")
    # not_last_rep = not_last_rep.drop(reseq_corrected_index.intersection(not_last_rep))  # type: ignore

    # drop the non terminal replicate
    unpaired_dataframe.loc[not_last_rep, "dropped_for"] = "non_terminal_replicate"
    remaining_dataframe = remaining_dataframe.drop(not_last_rep)
    logger.info(
        f"Step 0 - Dropping {len(not_last_rep):,} that are not last replicate, {len(last_rep_index):,} remaining"
    )

    # Step 1 - Drop samples with detected or expected controls
    expected_immo_index = remaining_dataframe[remaining_dataframe["expected_immo"]].index
    expected_synth_index = remaining_dataframe[remaining_dataframe["expected_synth"]].index
    expected_pooled_index = remaining_dataframe[remaining_dataframe["expected_pooled"]].index
    detected_immo_index = remaining_dataframe[remaining_dataframe["detected_immo"]].index
    detected_synth_index = remaining_dataframe[remaining_dataframe["detected_synth"]].index

    # log the indexes
    logger.info(f"Step 1 - Dropping {len(expected_immo_index):,} for failed expected immo   ")
    logger.info(f"Step 1 - Dropping {len(expected_synth_index):,} for failed expected synth")
    logger.info(f"Step 1 - Dropping {len(expected_pooled_index):,} for failed expected pooled")
    logger.info(f"Step 1 - Dropping {len(detected_immo_index):,} for failed detected immo   ")
    logger.info(f"Step 1 - Dropping {len(detected_synth_index):,} for failed detected synth")

    # mark why each one was dropped
    unpaired_dataframe.loc[detected_immo_index, "dropped_for"] = "detected_immo"
    unpaired_dataframe.loc[detected_synth_index, "dropped_for"] = "detected_synth"
    unpaired_dataframe.loc[expected_immo_index, "dropped_for"] = "expected_immo"
    unpaired_dataframe.loc[expected_synth_index, "dropped_for"] = "expected_synth"
    unpaired_dataframe.loc[expected_immo_index.intersection(expected_synth_index), "dropped_for"] = "expected_both"
    unpaired_dataframe.loc[expected_pooled_index, "dropped_for"] = "expected_pooled"

    # add all control indexes into one index seres
    control_indexes_to_drop = (
        expected_immo_index.append(expected_synth_index)  # type: ignore
        .append(expected_pooled_index)  # type: ignore
        .append(detected_immo_index)  # type: ignore
        .append(detected_synth_index)  # type: ignore
        .unique()  # type: ignore
    )  # type: ignore
    remaining_dataframe = remaining_dataframe.drop(control_indexes_to_drop)

    # Step 2 - Remove Empties which are negative controls or the mean sliced phred is na because we didn't detect a sequence
    expected_empty_index = remaining_dataframe[remaining_dataframe["is_neg_ctrl"]].index  # type: ignore
    detected_empty_index = remaining_dataframe[remaining_dataframe["mean_sliced_phred"].isna()].index  # type: ignore

    # drop detected empty before expected
    unpaired_dataframe.loc[detected_empty_index, "dropped_for"] = "detected_emtpy"
    unpaired_dataframe.loc[expected_empty_index, "dropped_for"] = "expected_empty"
    combined_empty_index = expected_empty_index.append(detected_empty_index).unique()  # type: ignore
    remaining_dataframe = remaining_dataframe.drop(combined_empty_index)  # type: ignore
    logger.info(f"Step 2 - Dropping {len(expected_empty_index):,} for expected emtpy")
    logger.info(
        f"Step 2 - Dropping {len(detected_empty_index.difference(expected_empty_index)):,} for detected empties"
    )

    # Step 3 - Remove Doublets
    expected_doublets_index = remaining_dataframe[remaining_dataframe["is_doublet"]].index
    logger.info(f"Step 3 - Dropping {len(expected_doublets_index):,} for expected doublets")
    unpaired_dataframe.loc[expected_doublets_index, "dropped_for"] = "expected_doublet"
    remaining_dataframe = remaining_dataframe.drop(expected_doublets_index)

    # Step 4 - Unproductive Mabs
    unproductive_index = remaining_dataframe[~remaining_dataframe["productive"].astype(bool)].index
    logger.info(f"Step 5 - Dropping {len(unproductive_index):,} for unproductive sequences with stop codons")
    unpaired_dataframe.loc[unproductive_index, "dropped_for"] = "unproductive_mab"
    remaining_dataframe = remaining_dataframe.drop(unproductive_index)

    # Step 5 - Incomplete Mabs
    incomplete_mabs_index = remaining_dataframe[~remaining_dataframe["complete_vdj"].astype(bool)].index
    logger.info(f"Step 4 - Dropping {len(incomplete_mabs_index):,} for incomplete antibodies")
    unpaired_dataframe.loc[incomplete_mabs_index, "dropped_for"] = "incomplete_mab"
    remaining_dataframe = remaining_dataframe.drop(incomplete_mabs_index)

    # Step 6 - Remove Bad Phred if its not corrected by reseq
    bad_phred_cutoff_index = remaining_dataframe.query(f"mean_sliced_phred < 30 & reseq_corrected==False").index
    logger.info(f"Step 6 - Dropping {len(bad_phred_cutoff_index):,} for phred below 30")
    unpaired_dataframe.loc[bad_phred_cutoff_index, "dropped_for"] = "low_phred"
    remaining_dataframe = remaining_dataframe.drop(bad_phred_cutoff_index)

    # Step 7 - Finally remove entries that can't be resolved by IGBLAST but aren't flagged a incomplete
    incomplete_by_software_mabs_index = remaining_dataframe[
        (
            (remaining_dataframe["fwr1_aa"].isna())
            | (remaining_dataframe["fwr2_aa"].isna())
            | (remaining_dataframe["fwr3_aa"].isna())
            | (remaining_dataframe["fwr4_aa"].isna())
            | (remaining_dataframe["cdr1_aa"].isna())
            | (remaining_dataframe["cdr2_aa"].isna())
            | (remaining_dataframe["cdr3_aa"].isna())
        )
    ].index
    logger.info(f"Step 7 - Dropping {len(incomplete_by_software_mabs_index):,} for incomplete by software lcdr3s")
    unpaired_dataframe.loc[incomplete_by_software_mabs_index, "dropped_for"] = "incomplete_by_software"
    remaining_dataframe = remaining_dataframe.drop(incomplete_by_software_mabs_index)
    # log the results so far
    dropped_for = unpaired_dataframe["dropped_for"].value_counts()  # type: ignore
    dropped_for["pairing_candidates"] = len(remaining_dataframe)
    logger.info(f"Pre-pairing results\n{dropped_for.sort_values()}")  # type: ignore

    # get all the things that don't have a mate...
    pairing_candidates = find_pairing_type(remaining_dataframe)
    bad_pair_index = pairing_candidates[~pairing_candidates["pairing_type"].isin(["IGH_IGK", "IGH_IGL"])].index

    # tell the unpaired dataframe about the bad pairs
    unpaired_dataframe.loc[bad_pair_index, "dropped_for"] = "unpaired"

    # the na is now that things that are IGH_IGK or IGH_IGL
    good_pairs_index = unpaired_dataframe[unpaired_dataframe["dropped_for"].isna()].index

    # get a series for logging purposes
    dropped_for_series = unpaired_dataframe["dropped_for"].value_counts()
    dropped_for_series["pairing_candidates"] = len(good_pairs_index)
    logger.info(f"Final results:\n{dropped_for_series}")
    return unpaired_dataframe


def pair(data: Data, working_dataframe: pd.DataFrame) -> Tuple[pd.DataFrame,pd.DataFrame]:
    """Take in the sequences and filter them to wells that had a single heavy light chain pair.
    This function willt take int he working dataframe and give back a paired dataframe with H/L pairs
    as well as a unpaired dataframe with the fates of each sequence.


    Parameters
    ----------
    data : Data
        The data object that contains the sequences and the metadata
    working_dataframe : pd.DataFrame
        The working dataframe from the NGS module
    Returns
    -------
    Tuple[pd.DataFrame]
        Two dataframes.
            paired_dataframe : The paired dataframe
            unpaired_dataframe : The unpaired dataframe with fates
    """
    logger.info(f"splitting {len(working_dataframe)} sequence_id data into model")
    unpaired_dataframe = drop_bad_samples(working_dataframe)
    paired_dataframe = pair_chains(unpaired_dataframe)
    return paired_dataframe, unpaired_dataframe
