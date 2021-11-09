#!/bin/bash
#PBS -l nodes=1:ppn=6:rhel7
#PBS -l walltime=09:00:00
#PBS -l pmem=8gb
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
wget -q https://datadryad.org/api/v2/datasets/doi%253A10.5061%252Fdryad.d5h62/download
unzip download
unzip nh_99-06.zip

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
cd Code
python 01_QC.py

echo ""
echo "-----PheEWAS analysis-----"
python 02_Analysis.py

echo ""
echo "-----Generate report-----"
jupyter nbconvert --execute --clear-output 03_ReportResults.ipynb 
jupyter nbconvert --to html 03_ReportResults.ipynb

echo "Job Ended at $(date)"
