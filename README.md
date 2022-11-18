# The IAVI G001 Clinical Trial Repository

[![Sequence Analysis Pipeline](https://github.com/SchiefLab/G001/workflows/Sequence%20Analysis%20Pipeline/badge.svg)](https://github.com/SchiefLab/G001/actions/workflows/integration.yml)
[![DOI](https://zenodo.org/badge/517925817.svg)](https://zenodo.org/badge/latestdoi/517925817)

- [The IAVI G001 Clinical Trial Repository](#the-iavi-g001-clinical-trial-repository)
- [Data Access](#data-access)
- [Pipeline](#pipeline)
  - [Installation Pre-requisites](#installation-pre-requisites)
  - [Installation](#installation)
  - [Flow Processing](#flow-processing)
  - [Collation of Flow Data](#collation-of-flow-data)
  - [Sequencing Pipeline](#sequencing-pipeline)
  - [Combine Flow and Sequence Data](#combine-flow-and-sequence-data)
  - [Figures](#figures)

# Data Access

If you don't want to run the code but would just like the important data files for your own analysis, you can find the following:

1. The raw flow cytometry files including .fcs, .xml and .csv files can be found [here](https://iavig001public.s3.us-west-2.amazonaws.com/flow_input.tgz)
   <br> **Warning - this file is large at ~120 GB and is hosted on a public Amazon S3 bucket**

2. The output of the processed are found in the [processed flow](data/flow/processed_flow/) directory.

3. The output of the processed flow needs to be collated from both trial sites into one [collated flow](data/flow/collated_flow/) directory.

4. The FASTQ files from Sanger sequencing are found in the [fastq](data/sequence/fastq) directory.

5. The annotated,filtered and paired antibody sequences are found in the [sequences/](data/figures/sequences/) directory and may be download with this [link](https://github.com/SchiefLab/G001/raw/main/data/figures/sequences/unblinded_sequences.csv.gz)

6. A merged summary file with all frequencies in this study can be found in the [flow_summary](data/figures/flow_summary) directory and may be download with this [link](https://github.com/SchiefLab/G001/raw/main/data/figures/flow_summary/flow_and_sequences.csv.gz)

# Pipeline

## Installation Pre-requisites

While not necessary, we highly recommend using the [conda](https://docs.conda.io/en/latest/) open-source package and environment manager. This allows you to make an environment with both Python and R dependencies. For the purposes of this repository, only minimal installer for anaconda is necessary (Miniconda).

<ins>Miniconda installers</ins>

[Mac command line installer](https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh)

[Mac GUI installer](https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.pkg)

[Linux command line installer](https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh)

Due to our dependencies on HMMER, there is no windows support at the moment.

## Installation

This installation assumes that `git` and `conda` are in your path.

```bash
# clone repository
git clone https://github.com/SchiefLab/G001.git

# change directory
cd G001

# install G001
./install.sh

# activate environment
conda activate G001

# initialize git-lfs
git-lfs install

# pull large files
git-lfs pull
```

## Flow Processing

The flow processing needs to be run on for each of the sites independently.

```bash
# get all raw flow files from public S3 bucket
wget https://iavig001public.s3.us-west-2.amazonaws.com/flow_input.tgz

# extract the files
tar -xvzf flow_input.tgz

# Run for FHCRC
g001 process-flow -s FHCRC -m FHCRC_manifest_file.csv -i flow_input/fhcrc/

# Run for VRC
g001 process-flow -s VRC -m VRC_manifest_file.csv -i flow_input/vrc/
```

## Collation of Flow Data

The following will combine the VRC and FHCRC flow data.

```bash
g001 collate -f data/flow/processed_flow/

# output will be in data/flow/processed_flow/Combined_Flow/
```

## Sequencing Pipeline

Analyze FASTQ sequences from Sanger sequencing.

```bash
G001 sequence_analysis -d data/ -r -o sequence_analysis_output
```

## Combine Flow and Sequence Data

```bash
Rscript src/g001/R/Combine_Flow_and_Seq_Results.R data/sequence data/flow/collated_flow
```

## Figures

The following can generate main text figures

```bash
g001 figures fig1
g001 figures fig2
g001 figures fig3
g001 figures fig4
g001 figures fig5
g001 figures fig6
g001 figures fig7
g001 figures fig8
```
