#!/bin/bash
#PBS -l walltime=5:00:00,nodes=1:ppn=1,mem=50gb
#PBS -o /home/hirschc1/della028/projects/sv_nams/analysis/qc
#PBS -e /home/hirschc1/della028/projects/sv_nams/analysis/qc
#PBS -V
#PBS -N summary_raw_gbs
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

module load R/3.6.0

# go to data folder
cd ~/projects/sv_nams/

# summarize data
Rscript scripts/summary_raw_gbs.R data/GBS-output/tmp analysis/qc/raw_gbs
