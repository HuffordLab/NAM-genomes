#!/bin/bash
#PBS -l walltime=12:00:00,nodes=1:ppn=1,mem=20gb
#PBS -o /home/hirschc1/della028/projects/sv_nams/data/GBS-output/tmp
#PBS -e /home/hirschc1/della028/projects/sv_nams/data/GBS-output/tmp
#PBS -V
#PBS -N filter_nam_rils_and_snps_within_svs_vcf_${CHR}
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

# go to project folder
cd ~/projects/sv_nams

{
  # skip header of "nam_ril_populations.txt" file
  read
  # set delimeter to tab
  IFS="\t"
  # read file line by line
  while read -r line; do
    # get name of the cross being parsed
    cross=$(echo $line | cut -f 1)
    # check if directory exists; if it doesnt, create one to store results
    [[ -d data/GBS-output/tmp/$cross ]] || mkdir -p data/GBS-output/tmp/$cross
    # transform line of the file into multiploe lines so that vcftools recognize 1 genotype to keep per line
    echo $line |  tr "\t" "\n" | tr "," "\n" > data/GBS-output/tmp/$cross/genotypes_to_keep.txt
    # use vcftools to filter a vcf file
    vcftools --vcf data/GBS-output/tmp/NAM_rils_SNPs.chr${CHR}.vcf \
             --keep data/GBS-output/tmp/$cross/genotypes_to_keep.txt \
             --exclude-bed data/tmp/SNPs_to_remove_$cross.bed \
             --out data/GBS-output/tmp/$cross/NAM_rils_SNPs.$cross.chr${CHR}.not-in-SVs \
             --recode \
             --recode-INFO-all
  done
} < "data/nam_ril_populations.txt"
