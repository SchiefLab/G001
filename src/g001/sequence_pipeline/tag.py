import numpy as np
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import pairwise_distances
from g001.data import Data
import pandas as pd
import logging
from pandarallel import pandarallel  # type: ignore
from Levenshtein._levenshtein import distance  # type: ignore

logger = logging.getLogger("tag")
np.random.seed(42)

def get_distance(x1: str, x2: str) -> int:
    x1 = str(x1)
    x2 = str(x2)
    return distance(x1, x2)  # type: ignore


def get_distance_norm(x1: str, x2: str) -> float:
    x1 = str(x1)
    x2 = str(x2)
    return 1 - (distance(x1, x2) / max(len(x1), len(x2)))  # type: ignore


def get_control_distances(sub_df_chains: pd.DataFrame, control_name: str, control_series: pd.Series) -> pd.DataFrame:  # type: ignore
    """
    get the distances between the control and the sub_df_chains
    """

    # we will check all these fields against the control fields
    check_fields = [
        "vdj_nt",
        "vdj_aa",
        "v_sequence_alignment",
        "cdr3",
        "cdr3_aa",
        "np1",
        "np2",
        "d_sequence_alignment",
        "j_sequence_alignment",
    ]

    def get_per_field_distance(fields: pd.Series, control_name: str, control_series: pd.Series) -> pd.Series:  # type: ignore
        all_series = []
        for control_metric in check_fields:
            distance = get_distance(fields[control_metric], control_series[control_metric])  # type: ignore
            distance_norm = get_distance_norm(fields[control_metric], control_series[control_metric])  # type: ignore
            all_series.append(  # type: ignore
                pd.Series(
                    {
                        f"{control_name.split('_')[0]}_{control_metric}_dist": distance,
                        f"{control_name.split('_')[0]}_{control_metric}_dist_norm": distance_norm,
                    }
                )
            )

        return pd.concat(all_series)  # type: ignore

    distance_df = (  # type: ignore
        sub_df_chains[check_fields]  # type: ignore
        .fillna("")  # type: ignore
        .parallel_apply(get_per_field_distance, axis=1, args=(control_name, control_series))
    )  # type: ignore
    return sub_df_chains.join(distance_df)  # type: ignore


def tag_controls(data: Data, working_dataframe: pd.DataFrame) -> pd.DataFrame:
    """Tag either Immo or Synth, the two conrols used in this experiment"""
    pandarallel.initialize(verbose=0)  # type: ignore

    # an airr representation of the conrol sequences used in G001
    controls_dataframe: pd.DataFrame = data.get_controls_datafrmae().set_index("sequence_id")  # type: ignore

    # first set that we didn't detect any controls
    working_dataframe.loc[:, "detected_synth"] = False
    working_dataframe.loc[:, "detected_immo"] = False

    # Now we only want to tag sequences with high mean quality scord
    working_dataframe_quality = working_dataframe[(working_dataframe["mean_sliced_phred"] >= 30)].copy()

    # parse out chains of working dataframe - these are detected chains
    heavy_chains = working_dataframe_quality[working_dataframe_quality["locus"] == "IGH"]
    kappa_chains = working_dataframe_quality[working_dataframe_quality["locus"] == "IGK"]
    lambda_chains = working_dataframe_quality[working_dataframe_quality["locus"] == "IGL"]

    # tag samples for IMMO hc,kc, lc
    logger.info(f"scanning {len(heavy_chains)} heavy chains for immo hc")
    sample_hc_vs_immo = get_control_distances(heavy_chains, "immo_hc", controls_dataframe.loc["immo_hc"])  # type: ignore
    logger.info(f"scanning {len(kappa_chains)} kappa chains for immo kc")
    sample_kc_vs_immo = get_control_distances(kappa_chains, "immo_kc", controls_dataframe.loc["immo_kc"])  # type: ignore
    logger.info(f"scanning {len(lambda_chains)} lambda chains for immo lc")
    sample_lc_vs_immo = get_control_distances(lambda_chains, "immo_lc", controls_dataframe.loc["immo_lc"])  # type: ignore

    # tag samples for SYNTH hc,kc, lc
    logger.info(f"scanning {len(heavy_chains)} heavy chains for synth hc")
    sample_hc_vs_synth = get_control_distances(heavy_chains, "synth_hc", controls_dataframe.loc["synth_hc"])  # type: ignore
    logger.info(f"scanning {len(kappa_chains)} kappa chains for immo kc")
    sample_kc_vs_synth = get_control_distances(kappa_chains, "synth_kc", controls_dataframe.loc["synth_kc"])  # type: ignore
    logger.info(f"scanning {len(lambda_chains)} lambda chains for synth lc")
    sample_lc_vs_synth = get_control_distances(lambda_chains, "synth_lc", controls_dataframe.loc["synth_lc"])  # type: ignore

    # Tagging criteria for samples
    # hc is 0.95 hcdr3 and 0.9 v gene
    immo_hc_tagging_index = sample_hc_vs_immo[
        (sample_hc_vs_immo["immo_cdr3_dist_norm"] > 0.95)
        & (sample_hc_vs_immo["immo_v_sequence_alignment_dist_norm"] > 0.9)
    ].index

    # norm distances of cdr3 of 0.92, np1 = 1 and v_seqence 0.9
    immo_kc_tagging_index = sample_kc_vs_immo[
        (sample_kc_vs_immo["immo_cdr3_dist_norm"] > 0.92)
        & (sample_kc_vs_immo["immo_v_sequence_alignment_dist_norm"] > 0.9)
        & (sample_kc_vs_immo["immo_np1_dist_norm"] == 1.0)
    ].index

    # norm distances of cdr3 of 0.92, and v_seqence 0.9
    immo_lc_tagging_index = sample_lc_vs_immo[
        (sample_lc_vs_immo["immo_cdr3_dist_norm"] > 0.92)
        & (sample_lc_vs_immo["immo_v_sequence_alignment_dist_norm"] > 0.9)
    ].index

    # tag immo with indexes for hc, kc and lc
    for taggers, tag_name in zip(
        [immo_hc_tagging_index, immo_kc_tagging_index, immo_lc_tagging_index], ["hc", "kc", "lc"]
    ):
        logger.info(f"tagging {len(taggers)} samples for immo {tag_name}")
        working_dataframe.loc[taggers, "detected_immo"] = True

    # Synth hc
    synth_hc_tagging_index = sample_hc_vs_synth[
        (sample_hc_vs_synth["synth_cdr3_dist_norm"] > 0.95)
        & (sample_hc_vs_synth["synth_v_sequence_alignment_dist_norm"] > 0.9)
    ].index

    # Synth kc
    synth_kc_tagging_index = sample_kc_vs_synth[
        (sample_kc_vs_synth["synth_cdr3_dist_norm"] > 0.9)
        & (sample_kc_vs_synth["synth_v_sequence_alignment_dist_norm"] > 0.9)
    ].index

    # Synth lc
    synth_lc_tagging_index = sample_lc_vs_synth[
        (sample_lc_vs_synth["synth_cdr3_dist_norm"] > 0.9)
        & (sample_lc_vs_synth["synth_v_sequence_alignment_dist_norm"] > 0.9)
    ].index

    # tag immo with indexes for hc, kc and lc
    for taggers, tag_name in zip(
        [synth_hc_tagging_index, synth_kc_tagging_index, synth_lc_tagging_index], ["hc", "kc", "lc"]
    ):
        logger.info(f"tagging {len(taggers)} samples for synth {tag_name}")
        working_dataframe.loc[taggers, "detected_synth"] = True
    return working_dataframe


def tag_intra_plate_duplicates(working_dataframe: pd.DataFrame, clust_distance: int) -> pd.DataFrame:
    """Cluster within each plate to find intra plate duplicates.

    Parameters
    ----------
    working_dataframe : pd.DataFrame
        working dataframe
    clust_distance : int
        The cutoff distance

    Returns
    -------
    pd.DataFrame
        working dataframe with intra plate duplicates tagged

    """
    logger.info("Setting up pandarallel to parallelize plate clustering")
    pandarallel.initialize(verbose=0, progress_bar=True)  # type: ignore

    group = ["pub_id", "plate", "timepoint", "chain", "replicate"]
    clustered_dataframe = working_dataframe.copy()
    lookup = ["vdj_nt"]

    def get_cluster_distance(df: pd.DataFrame) -> np.ndarray:  # type: ignore
        """Given a dataframe, get the N x N pairwise distances using Levenshtein distance of the lookup"""

        # need a lookup hash to get this quickly
        df_lookup = df[lookup].fillna("").to_dict(orient="index")  # type:ignore

        def calc_lev(x: np.array, y: np.array) -> int:  # type:ignore
            # calculate the Levenshtein distance between two strings
            string_1: str = df_lookup[x[0]]["vdj_nt"]  # type:ignore
            string_2: str = df_lookup[y[0]]["vdj_nt"]  # type:ignore
            if not isinstance(string_1, str):  # type:ignore
                string_1 = ""
            if not isinstance(string_2, str):  # type:ignore
                string_2 = ""
            if string_1 == "" or string_2 == "":  # type:ignore
                return clust_distance + 1  # will never let these be clustered
            return distance(string_1, string_2)  # type:ignore

        # turn that into a numpy array
        X = np.array(df.index).reshape(-1, 1)  # type:ignore
        return pairwise_distances(X, metric=calc_lev)  # type:ignore

    def apply_cluster(group_df: pd.DataFrame) -> pd.DataFrame:
        """
        The apply function that can be parallelized which will cluster each plate
        """

        # pandarallel is not finding this as a global variable so I'm redefining it here
        group = ["pub_id", "plate", "timepoint", "chain", "replicate"]

        # I have to take the cluster name out of the dataframe again
        group_name = list(group_df.groupby(group).groups.keys())  # type:ignore
        if len(group_name) > 1:
            raise ValueError("Duplicate group name")

        # now we can make the group name by on "_
        sub_group = "_".join([str(i) for i in group_name[0]])  # type:ignore

        # make a distance frame
        distance_frame = pd.DataFrame(
            get_cluster_distance(group_df),
            index=group_df.index,
            columns=group_df.index,
        )

        if len(distance_frame) == 1:
            # group_df.loc[:, "sub_cluster"] = "0"
            group_df.loc[:, "group_cluster"] = sub_group
            group_df.loc[:, "sub_cluster"] = 0
            group_df.loc[:, "intra_plate_total_cluster"] = sub_group + "_0"
            return group_df

        # create an agglorate frame
        model = AgglomerativeClustering(
            linkage="complete",
            n_clusters=None,
            affinity="precomputed",
            distance_threshold=clust_distance,
        )

        # cluster
        model.fit(distance_frame)  # type:ignore

        # get out name
        group_df.loc[:, "group_cluster"] = sub_group

        # get out the sub cluster within the plate
        group_df.loc[:, "sub_cluster"] = model.labels_  # type: ignore

        # join those two together to get a unique cluster name
        group_df.loc[:, "intra_plate_total_cluster"] = group_df[["group_cluster", "sub_cluster"]].apply(  # type:ignore
            lambda x: x[0] + "_" + str(x[1]), axis=1  # type:ignore
        )
        return group_df

    clustered_dataframe = clustered_dataframe.groupby(group).parallel_apply(apply_cluster)  # type:ignore
    return clustered_dataframe.reset_index(drop=True)  # type:ignore


def tag(data: Data, working_dataframe: pd.DataFrame) -> pd.DataFrame:
    """
    tag the working samples for controls and intra-plate duplicates
    """
    logger.info("Tagging intra-plate duplicates")
    working_dataframe = tag_intra_plate_duplicates(working_dataframe, 2).drop(["group_cluster", "sub_cluster"], axis=1)
    logger.info("Tagging possible controls")
    working_dataframe = tag_controls(data, working_dataframe)
    return working_dataframe
