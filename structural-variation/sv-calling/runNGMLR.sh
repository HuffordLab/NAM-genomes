#!/bin/bash
#genome="/work/LAS/mhufford-lab/arnstrm/newNAM/genomes/chromosomes-split-chr0/B73-chromosomes-split-chr0.fasta"
genome="$2"
query="$1"
out="$(basename ${genome%.*})_$(basename ${query%.*})"
echo $out
source /work/LAS/mhufford-lab/arnstrm/miniconda/etc/profile.d/conda.sh
conda activate ngmlr
ngmlr \
   --reference ${genome} \
   --query ${query} \
   --presets pacbio \
   --output ${out}.sam \
   --bam-fix \
   --no-progress \
   --threads 4

