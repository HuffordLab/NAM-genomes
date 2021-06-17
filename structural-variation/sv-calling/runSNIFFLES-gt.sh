#!/bin/bash
bamfile=$1
source /work/LAS/mhufford-lab/arnstrm/miniconda/etc/profile.d/conda.sh
conda activate ngmlr
ivcf=/work/LAS/mhufford-lab/arnstrm/newNAM/analyses/a-corrected-reads-mapping/SNIFFLES_v6_20200618/r1_vcf/survivor-merged_1kbpdist_r1.vcf
#sniffles \
#   --threads 12 \
#   --mapped_reads ${bamfile} \
#   --Ivcf ${ivcf} \
#   --bedpe ${bamfile%.*}_r2.bedpe &> ${bamfile%.*}_r2-sniffles-bed.log
sniffles \
   --threads 12 \
   --mapped_reads ${bamfile} \
   --Ivcf ${ivcf} \
   --vcf ${bamfile%.*}_r2.vcf &> ${bamfile%.*}_r2-sniffles-vcf.log
