#!/bin/bash
#PBS -l walltime=1:00:00,nodes=1:ppn=1,mem=100gb
#PBS -o /home/hirschc1/della028/projects/sv_nams/analysis/reseq_snps_projection2
#PBS -e /home/hirschc1/della028/projects/sv_nams/analysis/reseq_snps_projection2
#PBS -V
#PBS -N subset_hmp_based_on_ld_svs_${CHR}
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

# go to project folder
cd ~/projects/sv_nams/analysis/reseq_snps_projection2

run_pipeline.pl -Xmx120g -importGuess NAM_rils_SNPs-reseq_and_best-SNPs.reseq-snps-all-crosses.chr-${CHR}.projected.hmp.txt -includeSiteNamesInFile ~/projects/sv_nams/data/subset-NAM-snps/SVs-to-keep.missing-filter.no-tra.chr${CHR}.txt -export NAM_rils_subset_SVs.chr-${CHR}.hmp.txt -exportType HapmapDiploid
