from g001.data import Data
import pandas as pd
import logging

logger = logging.getLogger("unblind")


def unblind(data: Data, working_dataframe: pd.DataFrame) -> pd.DataFrame:
    """
    Add columns to the working dataframe such as the dose and placebo group.

    Also add the vrc01_class features
    """
    logger.info(f"unblinding {len(working_dataframe):,} pairs")

    # turn High Dose into high_dose, Low Dose into low_dose
    dose_group_df = data.get_dose_group().replace("Low Dose", "low_dose").replace("High Dose", "high_dose")

    # insert the vaccine group, placebo or vaccine
    working_dataframe.insert(
        6,
        "vaccine_group",
        working_dataframe["pub_id"].map(dose_group_df.set_index("pub_id")["vaccine_group"].to_dict()),
    )

    # insert the dose group, high or low
    working_dataframe.insert(
        6,
        "dose_group",
        working_dataframe["pub_id"].map(dose_group_df.set_index("pub_id")["dose_group"].to_dict()),
    )

    # does it have VH12 anywhere in the call?
    working_dataframe["has_vh1-2"] = (
        working_dataframe["v_call_heavy"]
        .fillna("")
        .str.split(",")
        .apply(lambda x: "IGHV1-2" in [i.split("*")[0] for i in x])
    )

    # insert cellid at the beginning that is the unique identifier
    working_dataframe.insert(
        0,
        "cellid",
        working_dataframe[["pub_id", "timepoint", "plate", "well"]].apply(
            lambda x: "_".join([i.replace("_", "") for i in x]), axis=1
        ),
    )

    # make sure cellid is unique
    if len(working_dataframe) != len(working_dataframe["cellid"].unique()):
        raise ValueError("cellid is not unique")

    # does it have 5 len
    working_dataframe["has_5_len_lcdr3"] = working_dataframe["cdr3_aa_light"].str.len() == 5

    # does it have both vh1-2 and 5 len
    working_dataframe["is_vrc01_class"] = working_dataframe["has_vh1-2"] & working_dataframe["has_5_len_lcdr3"]

    return working_dataframe
