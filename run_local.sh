#!/bin/bash

#### CREATE DIRECTORIES
mkdir Results
mkdir Results/Plots
mkdir Results/Plots/Inspect

#### RUN PIPELINE
echo '-----Downloading NHANES-----'
cd Data/
wget -q https://datadryad.org/stash/downloads/file_stream/70319
unzip 70319

cd ../Code/
echo '-----Running QC-----'
python 01_QC.py 
echo '-----Running Analysis-----'
python 02_Analysis.py
echo '-----Generate report-----'
jupyter nbconvert --execute --clear-output 03_ReportResults.ipynb 
jupyter nbconvert --to html 03_ReportResults.ipynb
