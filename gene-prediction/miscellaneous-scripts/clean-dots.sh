#!/bin/bash
ml bioawk
fasta=$1
bioawk -c fastx '{print $name}' $fasta > name.txt
bioawk -c fastx '{print $seq}' $fasta | sed 's/\.//g' > seq.txt
paste name.txt seq.txt | awk '{print ">"$1"\n"$2}' | fold > ${fasta}.1
rm name.txt seq.txt
