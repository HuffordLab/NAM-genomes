#!/bin/bash
#PBS -l walltime=12:00:00,nodes=1:ppn=5,mem=100gb
#PBS -o /home/hirschc1/della028/projects/sv_nams/analysis/reseq_snps_projection2
#PBS -e /home/hirschc1/della028/projects/sv_nams/analysis/reseq_snps_projection2
#PBS -V
#PBS -N add_mono_reseq_SNPs
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

# go to project folder
cd ~/projects/sv_nams/data/tmp/

module load parallel
module load R/3.6.0

parallel --jobs 5 < ~/projects/sv_nams/scripts/commands_add_mono-reseq-SNPs2.txt
