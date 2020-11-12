#!/bin/bash
source /work/LAS/mhufford-lab/arnstrm/miniconda/etc/profile.d/conda.sh
conda activate mikado
TransDecoder.LongOrfs \
   -t mikado_prepared.fasta
TransDecoder.Predict \
   -t mikado_prepared.fasta \
   --cpu 16

