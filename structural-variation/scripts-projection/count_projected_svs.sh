#!/bin/bash
#PBS -l walltime=8:00:00,nodes=1:ppn=10,mem=100gb
#PBS -o /home/hirschc1/della028/projects/sv_nams/analysis/reseq_snps_projection2
#PBS -e /home/hirschc1/della028/projects/sv_nams/analysis/reseq_snps_projection2
#PBS -V
#PBS -N count_projected_svs
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

# go to project folder
cd ~/projects/sv_nams/analysis/projection

module load R/3.6.0

Rscript ~/projects/sv_nams/scripts/count_projected_SVs.R ~/projects/sv_nams/data ~/projects/sv_nams/analysis/reseq_snps_projection2
