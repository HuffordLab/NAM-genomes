#!/bin/bash
#PBS -l walltime=1:00:00,nodes=1:ppn=1,mem=5gb
#PBS -o /home/hirschc1/della028/projects/sv_nams/data/GBS-output/tmp
#PBS -e /home/hirschc1/della028/projects/sv_nams/data/GBS-output/tmp
#PBS -V
#PBS -N vcf2hmp_nam-rils_snps-not-in-svs_${CROSS}
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

# go to project folder
cd ~/projects/sv_nams/data/GBS-output/tmp/

for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 scaffs; do
  vcf-sort ${CROSS}/NAM_rils_SNPs.${CROSS}.$chr.not-in-SVs.recode.vcf > ${CROSS}/NAM_rils_SNPs.${CROSS}.$chr.not-in-SVs.sorted.vcf
  run_pipeline.pl -Xmx10g -importGuess ${CROSS}/NAM_rils_SNPs.${CROSS}.$chr.not-in-SVs.sorted.vcf -export ${CROSS}/NAM_rils_SNPs.${CROSS}.$chr.not-in-SVs.hmp.txt -exportType HapmapDiploid
done
