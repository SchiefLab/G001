# G001

The G001 Clinical Trial Repository

[![Sequence Analysis Pipeline](https://github.com/SchiefLab/G001/workflows/Sequence%20Analysis%20Pipeline/badge.svg)](https://github.com/SchiefLab/G001/actions/workflows/integration.yml)
[![DOI](https://zenodo.org/badge/517925817.svg)](https://zenodo.org/badge/latestdoi/517925817)

# Downloading and Installing

You need git lfs in order to download the large files associated with this repository. You can download git lfs from [here](https://git-lfs.github.com/).

Once you have git lfs installed, you can clone the repository with the following command:

```bash
git clone https://github.com/SchiefLab/G001.git
cd G001
git lfs pull
```

To install

```bash
pip install .
```

# Sequencing Pipeline

```bash
G001 sequence_analysis -d data/
```
