from g001.data import Data
import pandas as pd
import logging

logger = logging.getLogger("Correct")


def correct(data: Data, working_dataframe: pd.DataFrame) -> pd.DataFrame:
    """
    Correct the data, unused..precalculated
    """
    logger.info(f"Corrected {len(working_dataframe):,}") 
    return working_dataframe
