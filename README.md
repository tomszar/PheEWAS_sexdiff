# PheEWAS sex differences

Analysis on sex differences in NHANES dataset. The analysis plan is [here](https://github.com/tomszar/PheEWAS_sexdiff/blob/main/ANALYSIS_PLAN.pdf).

## Clone repository

To clone the repository, on your terminal type:

```bash
git clone https://github.com/tomszar/PheEWAS_sexdiff.git
```

Then, enter the repository and follow the next instructions

```bash
cd PheEWAS_sexdiff
```

## Environment setup

The repository provides an `environment.yml` file to use with conda.

First, install [anaconda](https://www.anaconda.com/products/individual) on your local computer or user server account following the appropriate directions. Next, on the terminal, in the root of this repository, install the conda environment by running:

```bash
conda config --add channels r
conda config --add channels bioconda
conda config --add channels conda-forge
conda env create -f environment.yml
```

If the installation fails, you can install the environment manually, by running the following:

```bash
conda create --name py_clarite python=3.9 pandas numpy jupyterlab scipy statsmodels matplotlib scikit-learn r rpy2
conda activate py_clarite
pip install clarite
```

## Data

- The NHANES datasets were obtained from already unified dataset by [Patel et al 2016](https://doi.org/10/gdcc5d). The data is stored in [Dryad](https://datadryad.org/stash/dataset/doi:10.5061/dryad.d5h62)

## Replicate the analysis

To replicate the analysis, two bash scripts are provided; one is intended to run on a server, `run_server.sh`, while the other is intended to run in a local computer `run_local.sh`.
The `run_server.sh` file was ran at the Penn State [Roar system](https://www.icds.psu.edu/computing-services/roar-user-guide/), and used Roar's parameters.
If you are running the script on another server, you might want to tweak the parameters of the script.

- [Initial cleaning and QC](https://nbviewer.jupyter.org/github/tomszar/PheEWAS_sexdiff/blob/main/Code/02_QC.ipynb)
- [Sex difference analysis](https://nbviewer.jupyter.org/github/tomszar/PheEWAS_sexdiff/blob/main/Code/03_Analysis.ipynb)
