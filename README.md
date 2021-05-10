# PheEWAS sex differences

Analysis on sex differences in NHANES dataset. The analysis plan is [here](https://github.com/tomszar/PheEWAS_sexdiff/blob/main/ANALYSIS_PLAN.pdf).

## Data
- The NHANES datasets were obtained from already unified dataset by [Patel et al 2016](https://doi.org/10/gdcc5d). The data is stored in [Dryad](https://datadryad.org/stash/dataset/doi:10.5061/dryad.d5h62)
- [Script to download the data](https://github.com/tomszar/PheEWAS_sexdiff/blob/main/Code/01_GetData.sh)

## Requirements
- Create a conda environment with `conda env create -f environment.yml`. This step is already incorporated in the `run_repo.sh` script, so is not necessary.
- If there is a problem with the environment file, run this instead:
```
conda create --name py_clarite python=3.7 pandas numpy jupyterlab scipy statsmodels matplotlib
conda activate py_clarite
pip install git+https://github.com/HallLab/clarite-python.git
```

## Pipeline
- [Initial cleaning and QC](https://nbviewer.jupyter.org/github/tomszar/PheEWAS_sexdiff/blob/main/Code/02_QC.ipynb)
- [Sex difference analysis](https://nbviewer.jupyter.org/github/tomszar/PheEWAS_sexdiff/blob/main/Code/03_Analysis.ipynb): This was sent to Penn State's ICDS