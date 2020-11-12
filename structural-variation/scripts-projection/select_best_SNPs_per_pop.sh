#!/bin/bash
#PBS -l walltime=2:00:00,nodes=1:ppn=5,mem=50gb
#PBS -o /home/hirschc1/della028/projects/sv_nams/data/GBS-output/tmp
#PBS -e /home/hirschc1/della028/projects/sv_nams/data/GBS-output/tmp
#PBS -V
#PBS -N select_best_markers_per_pop_${CROSS}
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

module load R/3.6.0

# go to data folder
cd ~/projects/sv_nams/data/GBS-output/tmp/

# filter SNPs
Rscript ~/projects/sv_nams/scripts/select_best_SNPs_per_pop.R ${CROSS} ${CROSS}/NAM_rils_SNPs.${CROSS}.not-in-SVs.not-imputed.hmp.txt ${CROSS}/NAM_gbs-parents_SNPs.${CROSS}.not-in-SVs.reseq-overlay.hmp.txt ~/projects/sv_nams/analysis/qc/filter_best_SNPs --max_missing=0.3 --window_size=15 --window_step=1 --min_snps_per_window=5
