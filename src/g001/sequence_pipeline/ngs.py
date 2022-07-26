import logging
import warnings

import pandas as pd
from g001.data import Data
from Levenshtein._levenshtein import distance  # type: ignore
from sadie.airr import Airr

warnings.simplefilter("ignore", category=pd.errors.PerformanceWarning)  # type: ignore

logger = logging.getLogger("NGS")


def correct(ngs: pd.DataFrame, joined: pd.DataFrame) -> pd.DataFrame:
    """NGS correction algorithm"""

    # first, only get productive, complete from the NGS results
    ngs_result_dataframe_productive = (
        ngs[(ngs["productive"] & ngs["complete_vdj"] & ngs["vj_in_frame"])].copy().reset_index(drop=True)
    )

    # index we will use to join sanger and ngs dataframe
    index = ["pub_id", "timepoint", "plate", "chain", "well", "replicate"]

    # set an index to lookup
    joined_index = joined.set_index(index)  # type:ignore

    # this list wills store our ngs results that will update the sanger
    update_df = []

    # unfortuantely we have to go through each ngs index. There could be multilple contigs per index
    for row, group_df in ngs_result_dataframe_productive.groupby(  # type:ignore
        index
    ):
        try:
            single = joined_index.loc[row]
            if isinstance(single, pd.DataFrame):  # type:ignore
                raise ValueError(f"{row} in sanger dataframe has more than one entry")
        except KeyError:  # this happens when the index can't be found in ngs
            try:
                # try to find index again but this time, don't consider replicate
                single = joined_index.loc[row[:-1]]
            except KeyError:
                # we didn't find it even when dropping replicate
                logger.debug(f"{row[:-1]} missing everything in sanger")
                continue
            logger.debug(f"{row[-1]} missing in replicate in sanger but found {single.index[0]} in sanger for {row}")
            continue

        # if we found matching sanger and ngs....

        # first get the reference sequence_alignment of the vdj
        reference_seq_alignment: str = single["sequence_alignment"]

        # if it's None, we don't update because we just assumed nothing was in the sanger and not sure how it even got here
        if reference_seq_alignment is None or isinstance(reference_seq_alignment, float):
            logger.debug(f"skipping...{row} has no reference sequence_alignment")
            continue

        # if not, remove '-'
        reference_seq_alignment = reference_seq_alignment.replace("-", "")

        # for the ngs, go through all sequence alignment and find out how close it is to the sanger
        group_df["distance_to_ref"] = (
            group_df["sequence_alignment"]
            .str.replace("-", "")
            .apply(
                lambda x: distance(x, reference_seq_alignment)  # type:ignore
            )
        )

        # now find out if we have vdj_count > 100 and if th distnace is <=20 between reference sanger and ngs, we will correct it
        potential = group_df.query("vdj_count > 100 and distance_to_ref <=20")

        # if nothing matches, we yield to sanger
        if potential.empty:
            logger.debug(f"Nothing passed the filter, yielding to sanger {row}")
            continue
        else:
            # we get the sequence with the closest distance
            updatable = potential.sort_values("distance_to_ref").iloc[0]
            update_df.append(updatable)  # type:ignore

    # now we have collected all updatable and will turn it into a dataframe
    updatable = pd.DataFrame(update_df)

    # use the index to update the sanger with matching indexes from our updatable ngs
    joined_index.update(updatable.set_index(index))  # type:ignore

    # now log what we corrected
    joined_index["reseq_corrected"] = False

    # and set it to true for the ones that we corrected
    joined_index.loc[
        updatable.set_index(index).index,  # type:ignore
        "reseq_corrected",
    ] = True
    return joined_index.reset_index()  # reset index to get rid of multiindex


def ngs(data: Data, working_dataframe: pd.DataFrame) -> pd.DataFrame:
    """
    Annotate the NGS data and update ngs sequences where applicable
    """

    # these are the fields when combined are unique in the ngs data
    unique_grouping = ["pub_id", "timepoint", "plate", "chain", "well", "replicate", "rank"]
    ngs_dataframe = data.get_ngs_resequencing_file().drop(["vdj_sequence", "vdj_start", "vdj_end"], axis=1)

    # now lets make a unique name to have something to join back on later
    ngs_dataframe["unique_name"] = ngs_dataframe[unique_grouping].apply(  # type:ignore
        lambda x: "_".join([str(i) for i in x]), axis=1  # type:ignore
    )  # type: ignore

    # we want it to be unique so assert it is
    assert not ngs_dataframe["unique_name"].duplicated().any()  # type:ignore

    # use same annotation params as for the annotation module
    airr_api = Airr("human", adaptable=True)

    # run it as a dataframe on unique name and sequence
    ngs_results_dataframe = airr_api.run_dataframe(ngs_dataframe, "unique_name", "sequence")

    # merge it back on itself to get the original dataframe with annotations
    # also are changing sequence id to ngs_sequence id to distinguish from sanger seqeunce id.
    ngs_working_dataframe = (
        ngs_dataframe.drop("sequence", axis=1)
        .merge(  # type:ignore
            ngs_results_dataframe, left_on="unique_name", right_on="sequence_id"
        )
        .drop(["unique_name"], axis=1)
        .rename({"sequence_id": "ngs_sequence_id"}, axis=1)
    )

    # chain HEAVY/LIGHT to actual local
    ngs_working_dataframe["chain"] = ngs_working_dataframe["locus"].map(  # type:ignore
        {"IGH": "HEAVY", "IGK": "KAPPA", "IGL": "LAMBDA"}
    )

    # run the selection algorithm
    working_dataframe = correct(ngs_working_dataframe, working_dataframe)
    return working_dataframe
