<div align="center">
<img src="https://repository-images.githubusercontent.com/517925817/1eee6262-ea00-4269-82e3-c1950d5752d6" style="margin:0.5em;width:75%;height:75%">
</div>


[![Flow Process](https://github.com/SchiefLab/G001/actions/workflows/process.yml/badge.svg)](https://github.com/SchiefLab/G001/actions/workflows/process.yml)
[![Collate](https://github.com/SchiefLab/G001/actions/workflows/collate.yml/badge.svg)](https://github.com/SchiefLab/G001/actions/workflows/collate.yml)
[![Sequence Analysis Pipeline](https://github.com/SchiefLab/G001/workflows/Sequence%20Analysis%20Pipeline/badge.svg)](https://github.com/SchiefLab/G001/actions/workflows/integration.yml)
[![Combine](https://github.com/SchiefLab/G001/actions/workflows/combine.yml/badge.svg)](https://github.com/SchiefLab/G001/actions/workflows/combine.yml)
[![DOI](https://zenodo.org/badge/517925817.svg)](https://zenodo.org/badge/latestdoi/517925817)

This repository includes data and code used to produce the manuscript Leggat, Cohen, Willis, Fulp, deCamp et al. Science (2022). All data has been deidentified.


- [Data Access](#data-access)
- [Pipeline](#pipeline)
  - [Installation pre-requisites](#installation-pre-requisites)
  - [Installation](#installation)
  - [FACS analysis](#facs-analysis)
    - [Collation of flow data](#collation-of-flow-data)
  - [BCR sequence analysis](#bcr-sequence-analysis)
  - [Combined B cell frequency and BCR sequence analysis](#combined-b-cell-frequency-and-bcr-sequence-analysis)
- [Figures and tables](#figures-and-tables)
  - [Main figures](#main-figures)
  - [Tables](#tables)

# Data Access

If you don't want to run the code but would just like the important data files for your own analysis, you can find the following:

1. The raw FACS files including .fcs, .xml and .csv files can be found [here](https://iavig001public.s3.us-west-2.amazonaws.com/flow_input.tgz)
   <br> **Warning - this file is large at ~120 GB**

2. The outputs of the processed FACS files are found in the [processed flow](data/flow/flow_processed_out/) directory.

3. The outputs need to be collated from both trial sites into one [collated flow](data/flow/collated_flow/) directory.

4. The FASTQ files from Sanger sequencing are found in the [fastq](data/sequence/fastq) directory.

5. The annotated, filtered and paired antibody sequences are found in the [sequences/](data/figures/sequences/) directory and may be downloaded with this [link](https://github.com/SchiefLab/G001/raw/main/data/figures/sequences/unblinded_sequences.csv.gz)

6. A merged summary file with all frequencies reported in this study can be found in the [flow_summary](data/figures/flow_summary) directory and may be downloaded with this [link](https://github.com/SchiefLab/G001/raw/main/data/figures/flow_summary/flow_and_sequences.csv.gz)

# Pipeline

## Installation pre-requisites

While not necessary, we highly recommend using the [conda](https://docs.conda.io/en/latest/) open-source package and environment manager. This allows you to make an environment with both Python and R dependencies. For the purposes of this repository, only a minimal installer for anaconda is necessary (Miniconda).

<ins>Miniconda installers</ins>

[Mac command line installer](https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh)

[Mac GUI installer](https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.pkg)

[Linux command line installer](https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh)

Due to our dependencies on HMMER, there is no Windows support at the moment.

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
```

Optional - If you'd like to run the figure generation code, you must pull the large input files using `git-lfs`.

```bash
# initialize git-lfs
git-lfs install

# pull large files
git-lfs pull
```

Once you install you can use to list of available options
```bash
g001 --help
```

## FACS analysis

The flow processing needs to be run for the two sites (VRC,FHCRC) independently from the flow_input data. As of v2.0.0, this code is fully functional with deidentified data. 

```bash
# get all raw flow files from public S3 bucket
wget https://iavig001public.s3.us-west-2.amazonaws.com/flow_input.tgz

# extract the files
tar -xvzf flow_input.tgz

# Run for FHCRC
g001 process-flow -s fhcrc -i flow_input/fhcrc/ -o flow_processed_out/

# Run for VRC
g001 process-flow -s vrc -i flow_input/vrc/ -o flow_processed_out/

# For more options, use
g001 process-flow --help
```

### Collation of flow data

The following will combine the VRC and FHCRC flow data.

```bash
# If you ran the steps above in FACS analysis, you can use the following command to collate
g001 collate -i flow_processed_out/ -o collated_flow

# If you did not run the above steps in FACS analysis, we have run those steps and 
#placed the output in /data/flow/flow_process_out. You can use the following command to collate 
g001 collate -i data/flow/flow_process_out/ -o collated_flow

# For more options, use
g001 collate --help
```

## BCR sequence analysis

Run BCR sequence analysis pipeline (as in Leggat et al fig. S10). This code is fully functional with deidentified data.

```bash
# run sequence analysis and output to the folder sequence_analysis_output
g001 sequence_analysis -o sequence_analysis_output

# For more options, use
g001 sequence_analysis --help
```

## Combined B cell frequency and BCR sequence analysis

This code combines the sequencing and flow processing results and computes B cell frequencies among various sets of cells (e.g VRC01-class B cells among all IgG+ memory B cells). As of v2.0.0, this code is fully functional with deidentified data. 

```bash
# If you ran the steps above for collate and sequence analysis, 
# and the respective output folders are sequence_analysis_output and collated_flow, 
# you can use the following command to combine the sequence and flow data
g001 combine -s sequence_analysis_output -c collated_flow -o combined_flow_seq


# If you did not run the above steps for collate and sequence analysis, we have run those steps and placed the output 
# in /data/flow/collated_flow and data/sequence. You can use the following command to combine the sequence and flow data 
g001 combine -s data/sequence -c data/flow/collated_flow -o combined_flow_seq
```

# Figures and tables

## Main figures

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

## Tables

The following code generates all supplementary tables in the Leggat et al. manuscript. Both individual pdfs and a single combined pdf are generated. This command is only supported on Mac.

```
g001 supptables -c -o supp_tables
```





