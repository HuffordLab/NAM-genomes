#!/bin/bash
module load blast-plus
input=$1
blastx -max_target_seqs 5 -num_threads 4 -query ${input} -outfmt 5 -db uniprot-sprot_viridiplantae.fasta -evalue 0.000001 2> ${input%.*}.log | sed '/^$/d' | gzip -c - > ${input%.*}.blast.xml.gz
