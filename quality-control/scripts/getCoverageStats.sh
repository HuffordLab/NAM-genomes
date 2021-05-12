#!/bin/bash
file="$1"
nam=$(dirname $(dirname $file))
echo -e "NAM\tSample\tID\tTotalReads\tTotalAvgLen\tTotalMapped\tTotalMapped%\tCoverage\tEffCoverage"
grep -v "^Sample" $file | awk -v x=${nam} 'BEGIN{OFS=FS="\t"}{print x,$1,$2,$3,$4+$18,(($4+$18)/$2)*100,($2*$3)/169027888,($4*$3)/169027888}' |sed 's/_\(MN.....\)_R1/\t\1/g'
