#!/bin/bash
#PBS -l walltime=10:00:00,nodes=1:ppn=1,mem=110gb
#PBS -o /home/hirschc1/della028/projects/sv_nams/analysis/reseq_snps_projection2
#PBS -e /home/hirschc1/della028/projects/sv_nams/analysis/reseq_snps_projection2
#PBS -V
#PBS -N ld_distribution_snps-svs_${CHR}
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

# go to project folder
cd ~/projects/sv_nams/analysis/reseq_snps_projection2/

# transform hmp into plink format
run_pipeline.pl -Xmx110g -importGuess NAM_rils_SNPs-reseq_and_SVs-SNPs.reseq-snps-all-crosses.chr-${CHR}.projected.duplicated-SVs-removed.hmp.txt -export ld/NAM_rils_SNPs-reseq_and_SVs-SNPs.reseq-snps-all-crosses.chr-${CHR}.projected.duplicated-SVs-removed -exportType Plink

