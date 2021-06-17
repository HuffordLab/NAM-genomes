#!/bin/bash
bamfile=$1
source /work/LAS/mhufford-lab/arnstrm/miniconda/etc/profile.d/conda.sh
conda activate ngmlr
sniffles \
   --threads 12 \
   --mapped_reads ${bamfile} \
   --max_distance 1000 \
   --min_support 10 \
   --min_length 100 \
   --min_seq_size 5000 \
   --vcf ${bamfile%.*}_r1.vcf &> ${bamfile%.*}_snifflesb-vcf.log

