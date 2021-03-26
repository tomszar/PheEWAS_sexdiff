#!/bin/bash
#PBS -l nodes=1:ppn=1:rhel7
#PBS -l walltime=09:00:00
#PBS -l pmem=32gb
#PBS -A mah546_c_g_bc_default #or open
#PBS -j oe

#### THIS SCRIPT WILL INITIALIZE THE REPO BY:
# 1. CREATING THE RESULTS FOLDER
# 2. DOWNLOADING THE NHANES DATABASE
# 3. RUNNING THE ANALYSIS AND REPORTS

#### General notes
# This script will run for several hours depending on the system.
# You can run it in a local computer that you can leave on exclusively runnning this script, or you can send this script to a cluster that will run it
# This script used the Penn State roar server requirementThere is a .pbs file that you can take as an example to run the script in a cluster
# Consider that different systems have different requirements
# Finally, you can take this script as a suggestion and run it in different stepts

#### 1. CREATE FOLDER
echo "Job started on $(hostname) at $(date)"
#Go to working directory
cd ${PBS_O_WORKDIR}
echo ""
echo "------------ CREATING FOLDERS ------------"
mkdir Results
mkdir Results/Plots
mkdir Results/Plots/Inspect

#### 2. DOWNLOAD NHANES
echo ""
echo "------------ DOWNLOADING NHANES ------------"
cd Data
wget -q https://datadryad.org/stash/downloads/file_stream/71132
unzip 71132

#### 3. RUN ANALYSIS
echo ""
echo "------------ RUNNING ANALYSIS ------------"
echo ""
## 3.1. Activate conda environment
cd ..
echo "-----Setting up the conda environment-----"
echo ""
echo "Installing conda environment ... "
echo ""
conda env create -f environment.yml
echo "Activating conda environment ... "
echo ""
conda activate py_clarite

echo ""
echo "-----QC process-----"
jupyter nbconvert --execute --clear-output Code/02_QC.ipynb #Execute notebook
jupyter nbconvert --to html Code/02_QC.ipynb #Export it to html for visualization

echo ""
echo "-----PheEWAS analysis-----"
cd Code/
python 03_Analysis.py

echo ""
echo "-----Generate report-----"

echo "Job Ended at $(date)"
