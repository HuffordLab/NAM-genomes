#!/bin/bash
#PBS -l walltime=1:00:00,nodes=1:ppn=1,mem=50gb
#PBS -o /home/hirschc1/della028/projects/sv_nams/analysis/reseq_snps_projection2/
#PBS -e /home/hirschc1/della028/projects/sv_nams/analysis/reseq_snps_projection2/
#PBS -V
#PBS -N tassel_summary_projection_by_cross_${CROSS}
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

# go to project folder
cd ~/projects/sv_nams/analysis/reseq_snps_projection2

run_pipeline.pl -Xmx50g -importGuess ~/projects/sv_nams/analysis/reseq_snps_projection2/NAM_rils_SNPs-only.${CROSS}.poly.projected.hmp.txt -GenotypeSummaryPlugin -endPlugin -export ~/projects/sv_nams/analysis/reseq_snps_projection2/NAM_rils_SNPs-only_${CROSS}_OverallSummary
