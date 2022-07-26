import pandas as pd
from g001.data import Data
from g001.sequence_pipeline.anarci import run_mutational_analysis
from sadie.airr.airrtable import LinkedAirrTable


def mutation(data: Data, working_dataframe: pd.DataFrame) -> pd.DataFrame:
    """
    Add mutations_heavy and mutations_light
    """
    lat = LinkedAirrTable(working_dataframe, key_column="cellid")
    lat_with_mutational_analysis = run_mutational_analysis(lat, scheme="kabat")
    lat_with_mutational_analysis[["cellid", "mutations_heavy", "mutations_light"]]
    working_dataframe = working_dataframe.merge(
        lat_with_mutational_analysis[["cellid", "mutations_heavy", "mutations_light"]], on="cellid", how="outer"
    )
    return working_dataframe
