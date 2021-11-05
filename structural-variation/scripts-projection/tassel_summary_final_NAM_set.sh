#!/bin/bash
#PBS -l walltime=3:00:00,nodes=1:ppn=1,mem=100gb
#PBS -o /home/hirschc1/della028/projects/sv_nams/analysis/reseq_snps_projection2/ld
#PBS -e /home/hirschc1/della028/projects/sv_nams/analysis/reseq_snps_projection2/ld
#PBS -V
#PBS -N tassel_summary_${CHR}
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

# go to project folder
cd ~/projects/sv_nams/analysis/reseq_snps_projection2

run_pipeline.pl -Xmx100g -importGuess NAM_rils_SNPs-reseq_and_best-SNPs.reseq-snps-all-crosses.chr-${CHR}.projected.hmp.txt -GenotypeSummaryPlugin -endPlugin -export ld/missing_data_filter/tassel_summary_chr-${CHR}_NAM_rils_
