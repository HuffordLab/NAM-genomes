#!/bin/bash
source /work/LAS/mhufford-lab/arnstrm/miniconda/etc/profile.d/conda.sh
conda activate mikado
genome=$1
nam=$(basename $genome |cut -f 1 -d ".")
outdir=/work/LAS/mhufford-lab/arnstrm/newNAM/analyses/g-transcript-assembly/${nam}/mikado/portcullis_out
bam=/work/LAS/mhufford-lab/arnstrm/newNAM/analyses/g-transcript-assembly2/${nam}/merged_${nam}.bam

portcullis full \
   -t 36 \
   --use_csi \
   --output ${outdir} \
   --orientation FR \
   --strandedness firststrand \
     $genome \
     $bam
