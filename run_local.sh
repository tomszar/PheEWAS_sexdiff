#!/bin/bash

#### CREATE DIRECTORIES
mkdir Results
mkdir Results/Plots
mkdir Results/Plots/Inspect

#### RUN PIPELINE
echo '-----Downloading NHANES-----'
cd Data/
wget -q https://datadryad.org/api/v2/datasets/doi%253A10.5061%252Fdryad.d5h62/download
unzip download
unzip nh_99-06.zip

cd ../Code/
echo '-----Running QC-----'
python 01_QC.py 
echo '-----Running Analysis-----'
python 02_Analysis.py
echo '-----Generate report-----'
jupyter nbconvert --execute --clear-output 03_ReportResults.ipynb 
jupyter nbconvert --to html 03_ReportResults.ipynb
echo '-----Generate figures-----'
python 04_PlotFigures.py
