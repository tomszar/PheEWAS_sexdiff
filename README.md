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

To replicate the analysis, two bash scripts are provided; one is intended to run on a server, `run_server.sh`, while the other is intended to run on a local computer `run_local.sh`.

### On a server

On your server terminal, type

```bash
qsub run_server.sh
```

The parameters used in the script are the ones used in the Penn State Roar server. Depending on the system you will need to modify, remove, or add parameters to the script. The script also contains the creation and activation of the conda environment; you can alternatively make sure to install the environment before running the bash script and comment those lines.

### On a local computer

On the terminal type:

```bash
conda activate adni_project
bash run_local.sh > run_local.log
```
