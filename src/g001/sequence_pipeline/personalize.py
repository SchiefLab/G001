import logging

import pandas as pd
from g001.data import Data
from sadie.airr import Airr
from sadie.reference.reference import Reference, References

logger = logging.getLogger("personalize")


def personalize(data: Data, working_dataframe: pd.DataFrame) -> pd.DataFrame:
    """
    Split the data, unused..precalculated
    """
    haplotype_lookup = data.get_haplotype_data()

    # our baseline references is all available names
    baseline_references: References = References().from_yaml()

    # now we will get just human
    baseline_human: pd.DataFrame = baseline_references.get_dataframe().query("name=='human'")

    # take out VH12
    baseline_no_vh12 = baseline_human[baseline_human["gene"].str.split("*").str.get(0) != "IGHV1-2"].copy()

    references: References = References()
    haplotype_df = []

    # group by allele1 and 2
    for allele_group, allele_group_df in haplotype_lookup.groupby(["allele_1", "allele_2"]):

        # the haplotype name will be ex allele1_allele2, eg 02_02
        name = "_".join(list(map(lambda x: x.split("*")[-1], allele_group)))

        # only get the common
        need_alleles = list(set(list(allele_group)))

        # create a baseline ref with no vh12
        reference = Reference().from_dataframe(baseline_no_vh12.drop("name", axis=1))

        # add each allele one at a time
        for needed_allele in need_alleles:
            reference.add_gene({"species": "human", "gene": needed_allele, "source": "imgt"})
        logger.info(f"adding_references:{name}")

        # add the reference to references with name
        references.add_reference(name=name, reference=reference)
        allele_group_df["haplotype"] = name
        haplotype_df.append(allele_group_df)

    haplotype_df = pd.concat(haplotype_df).reset_index(drop=True)

    # now for each group, we need to extract the ptid subject
    personalized_df = []
    for haplo, haplot_df in haplotype_df.groupby("haplotype"):

        # only get ptids of that group
        pub_ids = haplot_df["pub_id"].to_list()

        # sub df is unblind dataframe with just the subjects of interest
        sub_df = working_dataframe[working_dataframe["pub_id"].isin(pub_ids)]

        # make an api call to get the allele
        airr_api = Airr(haplo, adaptable=True, references=references)
        heavy_airr_df = airr_api.run_dataframe(sub_df, "cellid", "sequence_heavy")
        light_airr_df = airr_api.run_dataframe(sub_df, "cellid", "sequence_light")

        # make some temporary dataframes to update
        h_ = heavy_airr_df.rename({"sequence_id": "cellid"}, axis=1).set_index("cellid")
        l_ = light_airr_df.rename({"sequence_id": "cellid"}, axis=1).set_index("cellid")

        # suffix with _heavy and light to update
        h_.columns = [i + "_heavy" for i in h_.columns]
        l_.columns = [i + "_light" for i in l_.columns]
        sub_df_cellid = sub_df.set_index("cellid")
        sub_df_cellid.update(h_)
        sub_df_cellid.update(l_)
        personalized_df.append(sub_df_cellid)

    # concat them all
    before_df_len = len(working_dataframe)
    working_dataframe = pd.concat(personalized_df).reset_index()
    if before_df_len != len(working_dataframe):
        raise ValueError(f"personalized {len(working_dataframe):,} != {before_df_len:,} before")
    return working_dataframe
