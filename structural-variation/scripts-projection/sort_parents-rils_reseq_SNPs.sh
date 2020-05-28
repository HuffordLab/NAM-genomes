#!/bin/bash
#PBS -l walltime=7:00:00,nodes=1:ppn=11,mem=55gb
#PBS -o /home/hirschc1/della028/projects/sv_nams/data/tmp
#PBS -e /home/hirschc1/della028/projects/sv_nams/data/tmp
#PBS -V
#PBS -N sort_reseq_SNPs
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

# go to project folder
cd ~/projects/sv_nams/data/tmp/

module load parallel
module load R/3.6.0

parallel --jobs 11 < ~/projects/sv_nams/scripts/commands_sort-parents_reseq-SNPs.txt
parallel --jobs 11 < ~/projects/sv_nams/scripts/commands_sort-rils_reseq-SNPs.txt
