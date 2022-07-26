from g001.data import Data
import pandas as pd
import logging

logger = logging.getLogger("join")


def get_sequence_id_of_manifest(row: pd.Series) -> str:  # type: ignore
    """
    get the sequence_id from the manifest row
    """
    pub_id: str = row["pub_id"].replace("_", "")  # type: ignore
    other = [str(i) for i in row[["timepoint", "plate", "chain", "well", "replicate"]]]  # type: ignore
    return "_".join([pub_id] + other)


def join(data: Data, working_dataframe: pd.DataFrame) -> pd.DataFrame:
    """
    join the sequences to the sequence manifest dataframe and return a combined dataframe
    """
    logger.info("Loading sequence manifest dataframes from VRC and FHCRC")
    vrc_manifest = data.get_vrc_flow_manifest()
    fhcrc_manifest = data.get_fhcrc_flow_manifest()
    logger.info("combining VRC and FHCRC dataframes")
    combined_manifest = pd.concat([vrc_manifest, fhcrc_manifest]).reset_index(drop=True)  # type: ignore
    combined_manifest["sequence_id"] = combined_manifest.apply(get_sequence_id_of_manifest, axis=1)  # type: ignore
    logger.info(
        f"merging sequence manifest {len(combined_manifest)} dataframe with working dataframe {len(working_dataframe)}"
    )
    joined_dataframe = combined_manifest.merge(working_dataframe.drop(["pub_id", "timepoint", "plate", "chain", "well", "replicate"], axis=1), on="sequence_id", how="inner")  # type: ignore
    logger.info(f"{len(joined_dataframe)} sequences joined")
    return joined_dataframe
