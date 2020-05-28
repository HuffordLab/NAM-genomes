#!/bin/bash

input=$1

diamond blastx \
   --threads 16 \
   --db uniref100 \
   --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovhsp\
   --out ${input%.*}-diamond.out \
   --query $input

