# The IAVI G001 Clinical Trial Repository

[![Sequence Analysis Pipeline](https://github.com/SchiefLab/G001/workflows/Sequence%20Analysis%20Pipeline/badge.svg)](https://github.com/SchiefLab/G001/actions/workflows/integration.yml)
[![DOI](https://zenodo.org/badge/517925817.svg)](https://zenodo.org/badge/latestdoi/517925817)

## Data Access

If you don't want to run the code but would just like the important data files for your own analysis, you can find the following: 

1. The raw flow cytometry files including .fcs, .xml and .csv files can be found [here](https://iavig001public.s3.us-west-2.amazonaws.com/raw_flow.tgz).
   <br> **Warning - this file is large and is thus hosted on a public Amazon S3 bucket**

2. The FASTQ files from Sanger sequencing are found in the [fastq](data/sequence/fastq) directory.

3. The annotated,filtered and paired antibody sequences are found in the [sequences/](data/figures/sequences/) directory and may be download with this [link](https://github.com/SchiefLab/G001/raw/main/data/figures/sequences/unblinded_sequences.csv.gz)

4. A merged summary file with all frequencies in this study can be found in the [flow_summary](data/figures/flow_summary) directory and may be download with this [link](https://github.com/SchiefLab/G001/raw/main/data/figures/flow_summary/flow_and_sequences.csv.gz)

## Downloading and Installing

<!-- You may need git lfs in order to download the large files associated with this repository. You can download git lfs from [here](https://git-lfs.github.com/). -->

<!-- Once you have git lfs installed, you can clone the repository with the following command: -->

```bash
# clone repository
git clone https://github.com/SchiefLab/G001.git

# change directory
cd G001

# create conda env
conda env create -f environment_cross_platform.yml --force
# activate the environment

conda activate G001
# install repository
pip install .

```
<!-- # git lfs pull -->


## Sequencing Pipeline

```bash
G001 sequence_analysis -d data/
```
