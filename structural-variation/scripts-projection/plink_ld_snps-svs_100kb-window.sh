#!/bin/bash
#PBS -l walltime=4:00:00,nodes=1:ppn=1,mem=100gb
#PBS -o /home/hirschc1/della028/projects/sv_nams/analysis/reseq_snps_projection2
#PBS -e /home/hirschc1/della028/projects/sv_nams/analysis/reseq_snps_projection2
#PBS -V
#PBS -N plink_ld_snps-svs_v5_${CHR}
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

# go to project folder
cd ~/projects/sv_nams/analysis/reseq_snps_projection2/

# calculate LD
plink --file ld/NAM_rils_SNPs-reseq_and_best-SNPs.reseq-snps-all-crosses.chr-${CHR}.projected.plk --make-founders --r2 gz --ld-window-r2 0 --ld-window 100000 --ld-window-kb 100 --out /scratch.global/della028/hirsch_lab/ld_files/NAM_rils_SNPs-reseq_and_best-SNPs.reseq-snps-all-crosses.chr-${CHR}.projected
