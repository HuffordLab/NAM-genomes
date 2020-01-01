#!/bin/bash
#source /work/LAS/mhufford-lab/arnstrm/programs/strawberry_1.1.1/sourceme
#module use /work/GIF/software/modules
#module load GIF/strawberry/1.1.1
PATH=$PATH:/work/LAS/mhufford-lab/arnstrm/programs/strawberry_1.1.1
bam="$1"
strawberry \
   --output-gtf strawberry-${bam%.*}.gtf \
   --logfile strawberry_assembled.log \
   --no-quant \
   --num-threads 36 \
   --verbose --fr \
   --min-transcript-size 100 \
     ${bam}
