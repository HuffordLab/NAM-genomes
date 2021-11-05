#!/bin/bash
#PBS -l walltime=4:00:00,nodes=1:ppn=1,mem=40gb
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

for chr in {1..10}; do
  Rscript ~/projects/sv_nams/scripts/merge_SVs_and_SNPs.R ${CROSS}/NAM_parents-reseq_SNPs.${CROSS}.chr${chr}.poly.not-in-SVs.hmp.txt ~/projects/sv_nams/data/NAM_parents_SVs-SNPs.${CROSS}.chr${chr}.sorted.hmp.txt ~/projects/sv_nams/data/NAM_rils_SVs-SNPs.${CROSS}.best-markers.not-projected.chr${chr}.sorted.hmp.txt NAM_parents_SNPs-reseq_and_SVs-SNPs.${CROSS}.poly.chr-${chr}.hmp.txt ${CROSS}/NAM_rils_SNPs-reseq_and_SVs-SNPs.${CROSS}.poly.chr-${chr}.not-projected.hmp.txt
done
