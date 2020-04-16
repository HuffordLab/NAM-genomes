#!/bin/bash
#PBS -l walltime=4:00:00,nodes=1:ppn=1,mem=80gb
#PBS -o /home/hirschc1/della028/projects/sv_nams/analysis/reseq_snps_projection2
#PBS -e /home/hirschc1/della028/projects/sv_nams/analysis/reseq_snps_projection2
#PBS -V
#PBS -N split_reseq-snps_by_chr_${CHR}
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

# go to project folder
cd ~/projects/sv_nams/analysis/reseq_snps_projection2

grep -v -P "^del|^dup|^ins|^inv|^tra" NAM_rils_SNPs-reseq_and_SVs-SNPs.reseq-snps-all-crosses.chr-${CHR}.projected.duplicated-SVs-removed.hmp.txt > NAM_rils_projected-reseq-SNPs-only.all-RILs.duplicated-SVs-removed.chr-${CHR}.hmp.txt

gzip -c NAM_rils_projected-reseq-SNPs-only.all-RILs.duplicated-SVs-removed.chr-${CHR}.hmp.txt > NAM_rils_projected-reseq-SNPs-only.all-RILs.duplicated-SVs-removed.chr-${CHR}.hmp.txt.gz
