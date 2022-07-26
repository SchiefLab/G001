from g001.data import Data
import pandas as pd
import logging
import warnings

logger = logging.getLogger("swap")


def swap(data: Data, working_dataframe: pd.DataFrame) -> pd.DataFrame:
    """
    correcting sample swap:

    PubID_164 : V08 and V10 is actually PubID_153: V08 and V10
    """
    logger.info(f"swapping PubID_164 -> PubID_153 at V08/V10 ")
    with warnings.catch_warnings():
        warnings.simplefilter(action="ignore", category=FutureWarning)
        update_1 = working_dataframe[
            (working_dataframe["pub_id"] == "PubID_164") & (working_dataframe["timepoint"].isin(["V08", "V10"]))  # type: ignore
        ].copy()
        update_1["pub_id"].replace({"PubID_164": "PubID_153"}, inplace=True)  # type: ignore
        logger.info(f"Swapping PubID_164 -> PubID_153 with {len(update_1)} samples")

        # Handle sample swap
        update_2 = working_dataframe[
            (working_dataframe["pub_id"] == "PubID_153") & (working_dataframe["timepoint"].isin(["V08", "V10"]))  # type: ignore
        ].copy()
        logger.info(f"Swapping PubID_153 -> PubID_164 with {len(update_2)} samples")
        update_2["pub_id"].replace({"PubID_153": "PubID_164"}, inplace=True)  # type: ignore

    # update in place
    working_dataframe.loc[:, "swapped"] = False
    working_dataframe.update(update_1, overwrite=True)  # type: ignore
    working_dataframe.loc[update_1.index, "swapped"] = True
    working_dataframe.update(update_2, overwrite=True)  # type: ignore
    working_dataframe.loc[update_2.index, "swapped"] = True
    return working_dataframe
