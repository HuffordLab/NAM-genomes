#!/bin/bash
source /work/LAS/mhufford-lab/arnstrm/miniconda/etc/profile.d/conda.sh
conda activate trinity
bam="$1"
out=$(basename ${bam%.*})
Trinity \
   --genome_guided_bam ${bam} \
   --max_memory 300G \
   --SS_lib_type FR \
   --min_contig_length 100 \
   --genome_guided_max_intron 10000 \
   --full_cleanup \
   --CPU 36
