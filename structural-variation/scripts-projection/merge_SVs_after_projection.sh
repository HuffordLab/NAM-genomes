#!/bin/bash
#PBS -l walltime=24:00:00,nodes=1:ppn=11,mem=100gb
#PBS -o /home/hirschc1/della028/projects/sv_nams/analysis/projection
#PBS -e /home/hirschc1/della028/projects/sv_nams/analysis/projection
#PBS -V
#PBS -N merge_SVs_after_projection
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

# go to project folder
cd ~/projects/sv_nams/analysis/reseq_snps_projection2

module load R/3.6.0

Rscript ~/projects/sv_nams/scripts/merge_SVs_after_projection.R ~/projects/sv_nams/data/NAM_founders_SVs.hmp.txt ~/projects/sv_nams/analysis/reseq_snps_projection2
