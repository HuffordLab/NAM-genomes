#!/bin/bash
#PBS -l walltime=3:00:00,nodes=1:ppn=1,mem=0gb
#PBS -o /home/hirschc1/della028/projects/sv_nams/analysis/reseq_snps_projection2
#PBS -e /home/hirschc1/della028/projects/sv_nams/analysis/reseq_snps_projection2
#PBS -V
#PBS -N plink_only_snp-sv-ld_${CHR}
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

# go to project folder
cd /scratch.global/della028/hirsch_lab/ld_files

# copy header for filtered file
zcat NAM_rils_SNPs-reseq_and_SVs-SNPs.reseq-snps-all-crosses.chr-${CHR}.projected.duplicated-SVs-removed_ld-w-100_v2.ld.gz | head -n 1 > NAM_rils_SNPs-reseq_and_SVs-SNPs.reseq-snps-all-crosses.chr-${CHR}.projected.duplicated-SVs-removed_ld-w-100_v2.no-tra.snp-sv.ld

# keep only snp and sv r2 (excluding translocations)
zcat NAM_rils_SNPs-reseq_and_SVs-SNPs.reseq-snps-all-crosses.chr-${CHR}.projected.duplicated-SVs-removed_ld-w-100_v2.ld.gz | awk '$3~/^del|^dup|^ins|^inv/ && $6~/^S/ || $3~/^S/ && $6~/^del|^dup|^ins|^inv/' - >> NAM_rils_SNPs-reseq_and_SVs-SNPs.reseq-snps-all-crosses.chr-${CHR}.projected.duplicated-SVs-removed_ld-w-100_v2.no-tra.snp-sv.ld
