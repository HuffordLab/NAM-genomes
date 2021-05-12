#!/bin/bash

counts="$1"
cut -f 1,7- $counts | grep -v "^#" | sed 's/_R1_round-2Aligned.sortedByCoord.out.bam//g' >  ${counts%.*}_cleaned.txt
awk 'NR>1' ${counts%.*}_cleaned.txt > .body
head -n 1 ${counts%.*}_cleaned.txt | tr "\t" "\n" |sed 's/\(^.*\)_\(MN.....\)/\2.\1/g' | tr "\n" "\t" | sed 's/\t$/\n/g' > .head
cat .head .body >> ${counts%.*}_cleaned-header.txt
base=$(dirname $counts)
sed 's/NAMLINE/'${base}'/g' DESeq2.R > ${base}/${base}_qc.R
chmod +x ${base}/${base}_qc.R
echo -e "id\tcondition" > ${base}/${base}_metadata
head -n 1 ${counts%.*}_cleaned.txt | tr "\t" "\n" |sed 's/\(^.*\)_\(MN.....\)/\2.\1\t\1/g' |grep -vi "geneid" >> ${base}/${base}_metadata
