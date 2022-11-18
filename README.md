# The IAVI G001 Clinical Trial Repository

[![Sequence Analysis Pipeline](https://github.com/SchiefLab/G001/workflows/Sequence%20Analysis%20Pipeline/badge.svg)](https://github.com/SchiefLab/G001/actions/workflows/integration.yml)
[![DOI](https://zenodo.org/badge/517925817.svg)](https://zenodo.org/badge/latestdoi/517925817)

This repository includes data and code used to produce the manuscript Leggat, Cohen, Willis, Fulp, deCamp et al. Science (2022)

- [The IAVI G001 Clinical Trial Repository](#the-iavi-g001-clinical-trial-repository)
- [Data Access](#data-access)
- [Pipeline](#pipeline)
  - [Installation pre-requisites](#installation-pre-requisites)
  - [Installation](#installation)
  - [FACS analysis](#facs-analysis)
    - [Collation of flow data](#collation-of-flow-data)
  - [BCR sequence analysis](#bcr-sequence-analysis)
  - [Combined B cell frequency and BCR sequence analysis](#combined-b-cell-frequency-and-bcr-sequence-analysis)
  - [Figures](#figures)

# Data Access

If you don't want to run the code but would just like the important data files for your own analysis, you can find the following:

1. The raw FACS files including .fcs, .xml and .csv files can be found [here](https://iavig001public.s3.us-west-2.amazonaws.com/flow_input.tgz)
   <br> **Warning - this file is large at ~120 GB**

2. The outputs of the processed FACS files are found in the [processed flow](data/flow/processed_flow/) directory.

3. The outputs need to be collated from both trial sites into one [collated flow](data/flow/collated_flow/) directory.

4. The FASTQ files from Sanger sequencing are found in the [fastq](data/sequence/fastq) directory.

5. The annotated,filtered and paired antibody sequences are found in the [sequences/](data/figures/sequences/) directory and may be downloaded with this [link](https://github.com/SchiefLab/G001/raw/main/data/figures/sequences/unblinded_sequences.csv.gz)

6. A merged summary file with all frequencies reported in this study can be found in the [flow_summary](data/figures/flow_summary) directory and may be downloaded with this [link](https://github.com/SchiefLab/G001/raw/main/data/figures/flow_summary/flow_and_sequences.csv.gz)

# Pipeline

## Installation pre-requisites

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

## FACS analysis

The flow processing needs to be run for the two sites (VRC,FHCRC) independently.

The code here was used to process the raw flow files during the trial. Those files contained information that needed to be deidentified before we could share the files themselves. The code as originally written was dependent on the format of the identifier data. Consequently, after we deidentified the raw data for inclusion in the repository, the original code was not able to process the deidentifed raw data files. We have included the code here for reference, but it is not functional on the deidentified data.

We are working to update the code to function with the deidentified data. Please check back for updates. In the meantime, we have included the deidentified versions of the original processed flow data used for trial analysis in the [processed flow](data/flow/processed_flow/) directory.

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

### Collation of flow data

The following will combine the VRC and FHCRC flow data.

```bash
g001 collate -f data/flow/processed_flow/

# output will be in data/flow/processed_flow/Combined_Flow/
```

## BCR sequence analysis

Run BCR sequence analysis pipeline (as in fig. S10)

```bash
g001 sequence_analysis -d data/ -r -o sequence_analysis_output
```

## Combined B cell frequency and BCR sequence analysis

This code combines the sequencing and flow processing results and computes B cell frequencies among various sets of cells (e.g VRC01-class B cells among all IgG+ memory B cells). Similar to the problem we experienced with [flow processing code](#facs-analysis), we found that the original code does not function properly on deidentified data.

As above, we are working to upgrade the code to work on deidentified data. Please check back for updates. In the meantime, we are providing deidentified versions of the output of the original code which can be found [here](data/figures/flow_summary/flow_and_sequences.csv.gz).

```bash
Rscript src/g001/R/Combine_Flow_and_Seq_Results.R data/sequence data/flow/collated_flow
```

## Figures

The following code generates the main text figures from the data in this repository.

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
