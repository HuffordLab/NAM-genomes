#!/bin/bash
module load stringtie
bam="$1"
out=$(basename ${bam%.*} )

stringtie \
   ${bam} \
   --rf \
   -m 100 \
   -p 36 \
   -v \
   -o ${out}_stringtie.gtf
