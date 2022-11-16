# The IAVI G001 Clinical Trial Repository

[![Sequence Analysis Pipeline](https://github.com/SchiefLab/G001/workflows/Sequence%20Analysis%20Pipeline/badge.svg)](https://github.com/SchiefLab/G001/actions/workflows/integration.yml)
[![DOI](https://zenodo.org/badge/517925817.svg)](https://zenodo.org/badge/latestdoi/517925817)

## Data Access

If you don't want to run the code but would just like the important data files for your own analysis, you can find the following:
1. The raw flow cytometry files including .fcs, .xml and .csv files can be found [here](https://iavig001public.s3.us-west-2.amazonaws.com/flow_input.tgz)
   <br> **Warning - this file is large at ~120 GB and is hosted on a public Amazon S3 bucket**

2. The output of the processed are found in the [processed flow](data/flow/processed_flow/) directory.

3. The FASTQ files from Sanger sequencing are found in the [fastq](data/sequence/fastq) directory.

4. The annotated,filtered and paired antibody sequences are found in the [sequences/](data/figures/sequences/) directory and may be download with this [link](https://github.com/SchiefLab/G001/raw/main/data/figures/sequences/unblinded_sequences.csv.gz)

5. A merged summary file with all frequencies in this study can be found in the [flow_summary](data/figures/flow_summary) directory and may be download with this [link](https://github.com/SchiefLab/G001/raw/main/data/figures/flow_summary/flow_and_sequences.csv.gz)

## Installation

While not necessary, we highly recommend using the conda open-source package and enviroment manager. This allows you to make an environment with both Python and R dependencies. For the purposes of this repostiory, only minimal installer for anaconda is necessary.

<ins>Miniconda installers</ins>

[Mac command line installer](https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh)

[Mac GUI installer](https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.pkg)

[Linux command line installer](https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh)

Due to our dependencies on HMMER, there is no windows support at the moment.

## Setup

This setup assumes that `git` and `conda` are in your path.

```bash
# clone repository
git clone https://github.com/SchiefLab/G001.git

# change directory
cd G001

# create conda env
conda env create -f environment_cross_platform.yml --force

# activate the environment
conda activate G001

# install G001
pip install .
```
## Flow Processing 

The flow processing needs to be run on for each of the sites independently.

```bash
# get all raw flow files from public S3 bucket
wget https://iavig001public.s3.us-west-2.amazonaws.com/flow_input.tgz

# extract the files
tar -xvzf flow_input.tgz

# Run for FHCRC
Rscript src/g001/R/Flow_Processing.R FHCRC flow_input/fhcrc/FHCRC_Flow_Manifest.csv flow_input/fhcrc/ yes flow_output

# Run for VRC
Rscript src/g001/R/Flow_Processing.R VRC flow_input/vrc/VRC_Flow_Manifest.csv flow_input/vrc/ yes flow_output
```

## Collation of Flow Data
The following will combine the VRC and FHCRC flow data.

```bash
Rscript src/g001/R/Collate_Flow_Data.R data/flow/processed_flow/
```


## Sequencing Pipeline

```bash
G001 sequence_analysis -d data/
```
