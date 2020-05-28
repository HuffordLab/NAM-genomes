#!/bin/bash
# fulltable
file1=$1
# can.1
file2=$2
awk 'BEGIN{OFS=FS="\t"} FNR==NR{a[$1]=$2 FS $0;next}{ print $0, a[$1]}' $file1 $file2 | \
sed 's/\t$//g' |\
sed 's/ /\t/g' |\
tr -s "\t" |\
cut -f 1-21,24- > ${file2%.*}_${file1%.*}.txt 
