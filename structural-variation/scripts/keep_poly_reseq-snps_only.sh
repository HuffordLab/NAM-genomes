#!/bin/bash
#PBS -l walltime=2:00:00,nodes=1:ppn=1,mem=40gb
#PBS -o /home/hirschc1/della028/projects/sv_nams/data/tmp
#PBS -e /home/hirschc1/della028/projects/sv_nams/data/tmp
#PBS -V
#PBS -N keep_poly_reseq-snps_only_${CROSS}
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

# go to project folder
cd ~/projects/sv_nams/data/tmp/

module load R

Rscript ~/projects/sv_nams/scripts/keep_poly_reseq-snps_only.R ${CROSS} ~/projects/sv_nams/data/tmp/${CROSS}/NAM_parents-reseq_SNPs.${CROSS}.not-in-SVs.hmp.txt
