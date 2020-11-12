#!/bin/bash
#PBS -l walltime=24:00:00,nodes=1:ppn=10,mem=130gb
#PBS -o /home/hirschc1/della028/projects/sv_nams/analysis/reseq_snps_projection2
#PBS -e /home/hirschc1/della028/projects/sv_nams/analysis/reseq_snps_projection2
#PBS -V
#PBS -N subsample_high-low_ld_${CHR}
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

# go to project folder
cd ~/projects/sv_nams/analysis/reseq_snps_projection2

module load R

Rscript ~/projects/sv_nams/scripts/subsample_SNPs_different_ld.R /scratch.global/della028/hirsch_lab/ld_files/NAM_rils_SNPs-reseq_and_best-SNPs.reseq-snps-all-crosses.chr-${CHR}.projected.no-tra.snp-sv.ld ~/projects/sv_nams/data/subset-NAM-snps/SVs-to-keep.missing-filter.no-tra.txt ~/projects/sv_nams/data/subset-NAM-snps/SNPs-to-keep.missing-filter.txt ld ~/projects/sv_nams/data/subset-NAM-snps ${CHR}
