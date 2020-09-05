#!/bin/bash
# merge gg markers and fasta file
if [ "$#" -ne 1 ] ; then
echo "please provide a NAM name for the marker file";
exit 0
fi
base=$1
cat *_chr*.fasta >> ${base}_merged.fasta
cat *.tsv >> ${base}_merged.txt

