#!/bin/bash
module load cufflinks
bam="$1"
out=$(basename ${bam%.*})
mkdir -p ${out}
cufflinks \
   --output-dir "${out}" \
   --num-threads 2 \
   --library-type fr-firststrand \
    --frag-len-mean 100 \
   --verbose \
   --no-update-check \
   ${bam}

