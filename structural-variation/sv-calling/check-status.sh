#!/bin/bash
nam=$1
if [ -f "${nam::-1}.done" ]; then
    echo -e "${nam::-1}\tDone\tDone\tDone"
else
cd $nam
current=$(ls *.sam | wc -l)
total=$(ls error_corrected_reads/*.fasta.gz | wc -l)
complete=$(cat *.log | awk '$7==0' | cut -f 9 | sort | wc -l)
echo -e "${nam::-1}\t$total\t$current\t$complete"
fi
