#!/bin/bash
#PBS -l walltime=24:00:00,nodes=1:ppn=24,mem=100gb
#PBS -o /home/hirschc1/della028/projects/sv_nams/analysis/reseq_snps_projection2
#PBS -e /home/hirschc1/della028/projects/sv_nams/analysis/reseq_snps_projection2
#PBS -V
#PBS -N sliding_window_reseq_SNPs_${CROSS}_${CHR}
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

# go to project folder
cd ~/projects/sv_nams/analysis/reseq_snps_projection2

module load R/3.6.0

Rscript ~/projects/sv_nams/scripts/sliding_window_approach_reseq-SNPs.R ${CROSS} ~/projects/sv_nams/analysis/reseq_snps_projection2/NAM_rils_SNPs-reseq_and_best-SNPs.${CROSS}.poly.chr-${CHR}.projected.hmp.txt ~/projects/sv_nams/analysis/reseq_snps_projection2/NAM_parents_SNPs-reseq_and_best-SNPs.${CROSS}.poly.chr-${CHR}.sorted.hmp.txt --window_size=45 --window_step=1 --min_snps_per_window=15
