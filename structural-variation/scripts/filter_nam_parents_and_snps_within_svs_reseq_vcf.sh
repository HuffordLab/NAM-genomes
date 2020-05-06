#!/bin/bash
#PBS -l walltime=30:00:00,nodes=1:ppn=1,mem=100gb
#PBS -o /home/hirschc1/della028/projects/sv_nams/data/tmp
#PBS -e /home/hirschc1/della028/projects/sv_nams/data/tmp
#PBS -V
#PBS -N filter_nam_parents_and_snps_within_svs_reseq_vcf
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
    echo "$cross"
    # check if directory exists; if it doesnt, create one to store results
    [[ -d data/tmp/$cross ]] || mkdir -p data/tmp/$cross
    # add the two parents of a cross in a file so that vcftools recognize 1 genotype to keep per line
    echo "B73" > data/tmp/$cross/genotypes_to_keep.txt
    echo $cross | cut -d "x" -f 2-3 >> data/tmp/$cross/genotypes_to_keep.txt
    # use vcftools to filter a vcf file
    vcftools --vcf data/tmp/NAM_founders_SNPs.vcf \
             --keep data/tmp/$cross/genotypes_to_keep.txt \
             --exclude-bed data/tmp/SNPs_to_remove_$cross.bed \
             --out data/tmp/$cross/NAM_parents-reseq_SNPs.$cross.not-in-SVs \
             --recode
  done
} < "data/nam_ril_populations.txt"
