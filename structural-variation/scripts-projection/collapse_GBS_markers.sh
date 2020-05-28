#!/bin/bash
#PBS -l walltime=20:00:00,nodes=1:ppn=10,mem=120gb
#PBS -o /home/hirschc1/della028/projects/sv_nams/data/GBS-output/tmp
#PBS -e /home/hirschc1/della028/projects/sv_nams/data/GBS-output/tmp
#PBS -V
#PBS -N collapse_GBS_markers
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

module load R/3.6.0
module load parallel

# go to data folder
cd ~/projects/sv_nams/data/GBS-output/tmp/

# filter SNPs
parallel --tmpdir /scratch.global/della028/hirsch_lab/sv_nams/data/GBS-output/tmp/ --jobs 10 < ~/projects/sv_nams/scripts/commands_for_collapse-GBS-SNPs.txt
