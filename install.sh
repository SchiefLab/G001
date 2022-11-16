# determins bash, zsh or fish
eval "$(conda shell.bash hook)"
conda env create -f environment_cross_platform.yml --force
# conda create -n G001 
conda activate G001
# conda install -y r-remotes

Rscript -e "remotes::install_github('thebioengineer/colortable', dependencies = FALSE)"
Rscript -e "remotes::install_github('RGLab/cytoqc', dependencies = FALSE)"

pip install .