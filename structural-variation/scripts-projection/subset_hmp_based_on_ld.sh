#!/bin/bash
#PBS -l walltime=2:00:00,nodes=1:ppn=1,mem=120gb
#PBS -o /home/hirschc1/della028/projects/sv_nams/analysis/reseq_snps_projection2
#PBS -e /home/hirschc1/della028/projects/sv_nams/analysis/reseq_snps_projection2
#PBS -V
#PBS -N subset_hmp_based_on_ld_${CHR}
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

# go to project folder
cd ~/projects/sv_nams/analysis/reseq_snps_projection2

run_pipeline.pl -Xmx120g -importGuess NAM_rils_SNPs-reseq_and_best-SNPs.reseq-snps-all-crosses.chr-${CHR}.projected.hmp.txt -includeSiteNamesInFile ld/subset_low-ld-snps/SNPs-kept_chr${CHR}.txt -export NAM_rils_subset_SNPs.chr-${CHR}.snps-low-ld-sv.hmp.txt -exportType HapmapDiploid

run_pipeline.pl -Xmx120g -importGuess NAM_rils_SNPs-reseq_and_best-SNPs.reseq-snps-all-crosses.chr-${CHR}.projected.hmp.txt -includeSiteNamesInFile ld/subset_high-ld-snps/SNPs-kept_chr${CHR}.txt -export NAM_rils_subset_SNPs.chr-${CHR}.snps-high-ld-sv.hmp.txt -exportType HapmapDiploid

run_pipeline.pl -Xmx120g -importGuess NAM_rils_SNPs-reseq_and_best-SNPs.reseq-snps-all-crosses.chr-${CHR}.projected.hmp.txt -includeSiteNamesInFile ld/subset_random-snps/SNPs-kept_chr${CHR}.txt -export NAM_rils_subset_SNPs.chr-${CHR}.snps-random.hmp.txt -exportType HapmapDiploid
