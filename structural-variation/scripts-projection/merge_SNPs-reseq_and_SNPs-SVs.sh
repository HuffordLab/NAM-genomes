#!/bin/bash
#PBS -l walltime=4:00:00,nodes=1:ppn=1,mem=60gb
#PBS -o /home/hirschc1/della028/projects/sv_nams/data/tmp
#PBS -e /home/hirschc1/della028/projects/sv_nams/data/tmp
#PBS -V
#PBS -N merge_SNPs-reseq_and_SNPs-SVs_${CROSS}
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

# go to project folder
cd ~/projects/sv_nams/data/tmp/

module load R

Rscript ~/projects/sv_nams/scripts/merge_SNPs-reseq_and_SNPs-SVs.R ~/projects/sv_nams/data/tmp/${CROSS}/NAM_parents-reseq_SNPs.${CROSS}.poly.not-in-SVs.hmp.txt ~/projects/sv_nams/data/NAM_parents_SVs-SNPs.${CROSS}.sorted.hmp.txt ~/projects/sv_nams/analysis/projection/NAM_rils_SVs-SNPs.${CROSS}.best-markers.projected.hmp.txt ~/projects/sv_nams/data/tmp/NAM_parents_SNPs-reseq_and_SVs-SNPs.${CROSS}.poly.hmp.txt ~/projects/sv_nams/data/tmp/${CROSS}/NAM_rils_SNPs-reseq_and_SVs-SNPs.${CROSS}.poly.not-projected.hmp.txt
