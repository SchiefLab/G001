# determins bash, zsh or fish
eval "$(conda shell.bash hook)"
conda env create -f environment_cross_platform.yml --force
conda activate G001

# R install from github directly (not avaliable conda)
Rscript -e "remotes::install_github('thebioengineer/colortable', dependencies = FALSE)"
Rscript -e "remotes::install_github('RGLab/cytoqc', dependencies = FALSE)"

pip install -e .