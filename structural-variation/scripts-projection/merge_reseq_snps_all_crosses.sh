#!/bin/bash
#PBS -l walltime=8:00:00,nodes=1:ppn=1,mem=100gb
#PBS -o /home/hirschc1/della028/projects/sv_nams/data/tmp
#PBS -e /home/hirschc1/della028/projects/sv_nams/data/tmp
#PBS -V
#PBS -N merge_reseq_snps_all_crosses
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

# go to project folder
cd ~/projects/sv_nams/data/tmp/

# merge all crosses with projected resequencing SNPs
for chr in {1..10}; do
  # exclude first columns for all crosses
  cd ~/projects/sv_nams/data/tmp/
  for cross in $(ls -d B73x*); do
    echo $chr $cross
    cut -f 1-11 --complement ~/projects/sv_nams/analysis/reseq_snps_projection2/NAM_rils_SNPs-only.$cross.reseq-snps-all-crosses.chr-$chr.projected.hmp.txt > ~/projects/sv_nams/analysis/reseq_snps_projection2/NAM_rils_SNPs-only.$cross.reseq-snps-all-crosses.chr-$chr.rils-only.projected.hmp.txt
  done
  # join all rils in one file (keep entire hmp file for cross B73xB97 though)
  cd ~/projects/sv_nams/analysis/reseq_snps_projection2
  paste NAM_rils_SNPs-only.B73xB97.reseq-snps-all-crosses.chr-$chr.projected.hmp.txt \
        NAM_rils_SNPs-only.B73xCML103.reseq-snps-all-crosses.chr-$chr.rils-only.projected.hmp.txt \
        NAM_rils_SNPs-only.B73xCML228.reseq-snps-all-crosses.chr-$chr.rils-only.projected.hmp.txt \
        NAM_rils_SNPs-only.B73xCML247.reseq-snps-all-crosses.chr-$chr.rils-only.projected.hmp.txt \
        NAM_rils_SNPs-only.B73xCML277.reseq-snps-all-crosses.chr-$chr.rils-only.projected.hmp.txt \
        NAM_rils_SNPs-only.B73xCML322.reseq-snps-all-crosses.chr-$chr.rils-only.projected.hmp.txt \
        NAM_rils_SNPs-only.B73xCML333.reseq-snps-all-crosses.chr-$chr.rils-only.projected.hmp.txt \
        NAM_rils_SNPs-only.B73xCML52.reseq-snps-all-crosses.chr-$chr.rils-only.projected.hmp.txt \
        NAM_rils_SNPs-only.B73xCML69.reseq-snps-all-crosses.chr-$chr.rils-only.projected.hmp.txt \
        NAM_rils_SNPs-only.B73xHp301.reseq-snps-all-crosses.chr-$chr.rils-only.projected.hmp.txt \
        NAM_rils_SNPs-only.B73xIl14H.reseq-snps-all-crosses.chr-$chr.rils-only.projected.hmp.txt \
        NAM_rils_SNPs-only.B73xKi11.reseq-snps-all-crosses.chr-$chr.rils-only.projected.hmp.txt \
        NAM_rils_SNPs-only.B73xKi3.reseq-snps-all-crosses.chr-$chr.rils-only.projected.hmp.txt \
        NAM_rils_SNPs-only.B73xKy21.reseq-snps-all-crosses.chr-$chr.rils-only.projected.hmp.txt \
        NAM_rils_SNPs-only.B73xM162W.reseq-snps-all-crosses.chr-$chr.rils-only.projected.hmp.txt \
        NAM_rils_SNPs-only.B73xM37W.reseq-snps-all-crosses.chr-$chr.rils-only.projected.hmp.txt \
        NAM_rils_SNPs-only.B73xMo18W.reseq-snps-all-crosses.chr-$chr.rils-only.projected.hmp.txt \
        NAM_rils_SNPs-only.B73xMS71.reseq-snps-all-crosses.chr-$chr.rils-only.projected.hmp.txt \
        NAM_rils_SNPs-only.B73xNC350.reseq-snps-all-crosses.chr-$chr.rils-only.projected.hmp.txt \
        NAM_rils_SNPs-only.B73xNC358.reseq-snps-all-crosses.chr-$chr.rils-only.projected.hmp.txt \
        NAM_rils_SNPs-only.B73xOh43.reseq-snps-all-crosses.chr-$chr.rils-only.projected.hmp.txt \
        NAM_rils_SNPs-only.B73xOh7B.reseq-snps-all-crosses.chr-$chr.rils-only.projected.hmp.txt \
        NAM_rils_SNPs-only.B73xP39.reseq-snps-all-crosses.chr-$chr.rils-only.projected.hmp.txt \
        NAM_rils_SNPs-only.B73xTx303.reseq-snps-all-crosses.chr-$chr.rils-only.projected.hmp.txt \
        NAM_rils_SNPs-only.B73xTzi8.reseq-snps-all-crosses.chr-$chr.rils-only.projected.hmp.txt > NAM_rils_projected-SNPs-SVs.reseq-snps-all-crosses.chr-$chr.projected.hmp.txt
done
