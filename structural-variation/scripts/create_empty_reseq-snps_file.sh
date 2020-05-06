#!/bin/bash
#PBS -l walltime=3:00:00,nodes=1:ppn=1,mem=40gb
#PBS -o /home/hirschc1/della028/projects/sv_nams/analysis/reseq_snps_projection2
#PBS -e /home/hirschc1/della028/projects/sv_nams/analysis/reseq_snps_projection2
#PBS -V
#PBS -N create_empty_reseq_snps
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

# go to project folder
cd ~/projects/sv_nams/analysis/reseq_snps_projection2

module load R/3.6.0

Rscript ~/projects/sv_nams/scripts/create_empty_all-reseq-snps_file.R ${CHR} ~/projects/sv_nams/analysis/reseq_snps_projection2
