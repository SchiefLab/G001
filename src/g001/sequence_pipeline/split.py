from g001.data import Data
import pandas as pd
import logging

logger = logging.getLogger("Split")


def split(data: Data, working_dataframe: pd.DataFrame) -> pd.DataFrame:
    """
    Split the data, unused..precalculated
    """
    logger.info(f"splitting {len(working_dataframe):,} sequence_id data into model")
    return working_dataframe
