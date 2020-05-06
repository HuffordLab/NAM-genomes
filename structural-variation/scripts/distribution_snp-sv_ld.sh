#!/bin/bash
#PBS -l walltime=4:00:00,nodes=1:ppn=10,mem=120gb
#PBS -o /home/hirschc1/della028/projects/sv_nams/analysis/reseq_snps_projection2
#PBS -e /home/hirschc1/della028/projects/sv_nams/analysis/reseq_snps_projection2
#PBS -V
#PBS -N distribution_snp-sv_ld_${CHR}
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

# go to project folder
cd ~/projects/sv_nams/analysis/reseq_snps_projection2

module load R

# low ld
Rscript ~/projects/sv_nams/scripts/distribution_snps-LD-svs.R ld/plink_results_SNPs-lowest-LD-SV_chr${CHR}.missing-filter.ld ld/subset_low-ld-snps ~/projects/sv_nams/data/subset-NAM-snps/SNPs-to-keep_subsample-low-ld_chr${CHR}.missing-filter.txt

# high ld
Rscript ~/projects/sv_nams/scripts/distribution_snps-LD-svs.R ld/plink_results_SNPs-highest-LD-SV_chr${CHR}.missing-filter.ld ld/subset_high-ld-snps ~/projects/sv_nams/data/subset-NAM-snps/SNPs-to-keep_subsample-high-ld_chr${CHR}.missing-filter.txt

# random snps
Rscript ~/projects/sv_nams/scripts/distribution_snps-LD-svs.R /scratch.global/della028/hirsch_lab/ld_files/NAM_rils_SNPs-reseq_and_best-SNPs.reseq-snps-all-crosses.chr-${CHR}.projected.no-tra.snp-sv.ld ld/subset_random-snps ~/projects/sv_nams/data/subset-NAM-snps/SNPs-to-keep_random-missing-filter_chr${CHR}.txt
