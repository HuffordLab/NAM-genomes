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

run_pipeline.pl -Xmx100g -importGuess NAM_rils_projected-SNPs-only.all-RILs.chr-${CHR}.v9.hmp.txt.gz -GenotypeSummaryPlugin -endPlugin -export tassel_summary_NAM_rils_chr-${CHR}_projected_snps
