#!/bin/bash
for fq in $(find . -name "*_R1.fq.gz" | sed 's/\.\///g'); do
echo $fq |sed 's/_MN....._R1.fq.gz//g' |sed 's/\//\t/g' ;
done | datamash --sort crosstab 1,2 |\
sed 's/N\/A/0/g'
