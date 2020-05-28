#!/bin/bash
#PBS -l walltime=4:00:00,nodes=1:ppn=1,mem=80gb
#PBS -o /home/hirschc1/della028/projects/sv_nams/analysis/reseq_snps_projection2
#PBS -e /home/hirschc1/della028/projects/sv_nams/analysis/reseq_snps_projection2
#PBS -V
#PBS -N sort_reseq_snps_and_svs_chr_${CHR}
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

# go to project folder
cd ~/projects/sv_nams/analysis/reseq_snps_projection2

run_pipeline.pl -Xmx100g -SortGenotypeFilePlugin -inputFile NAM_rils_SNPs-reseq_and_best-SNPs.reseq-snps-all-crosses.chr-${CHR}.projected.hmp.txt -outputFile NAM_rils_SNPs-reseq_and_best-SNPs.reseq-snps-all-crosses.chr-${CHR}.projected.hmp.txt -fileType Hapmap
