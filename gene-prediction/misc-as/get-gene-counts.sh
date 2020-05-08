#!/bin/bash
res="phylostrata_table.csv"
out=$(echo $res | sed 's/_table.csv/_gene-stats.tsv/g')
cut -f 1,2 -d "," phylostrata_stats.csv |sed -e 's/"//g' -e 's/,/\t/g' -e 's/ /_/g' | grep -v "Var1" | awk '{print $2"\t"$1}' > .order
sed -e 's/,/\t/g' -e 's/"//g' -e 's/ /_/g' $res |\
cut -f 2,5 |\
sed 's/_P/\tP/1' |\
cut -f 1,3 |\
sort -k1,1 -u |\
cut -f 2 |grep -v "qseqid" | sort | uniq -c | awk '{print $2"\t"$1}' > .counts
awk 'BEGIN{OFS=FS="\t"} FNR==NR{a[$1]=$2 FS $3;next}{ print $0, a[$1]}' .order .counts |sort -k3,3 -n |cut -f 1,2 > $out
rm .counts .order

