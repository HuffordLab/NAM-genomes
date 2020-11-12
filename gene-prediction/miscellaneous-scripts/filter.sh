#!/bin/bash
aed=$(find $(pwd) -name "*gff.aed.filtered")
ids=$(find $(pwd) -name "*_primary-ids.txt")
res=phylostrata_table.csv
out=gene-stats_filtered.tsv
out2=gene-stats.tsv
sed 's/_T/\tT/1' $aed | cut -f 1 |sort |uniq > .ids
cut -f 1,2 -d "," phylostrata_stats.csv | sed -e 's/"//g' -e 's/,/\t/g' -e 's/ /_/g' | grep -v "Var1" | awk '{print $2"\t"$1}' > .order
sed -e 's/,/\t/g' -e 's/"//g' -e 's/ /_/g' $res |\
cut -f 2,5 |\
sed 's/_P/_T/1' |\
grep -v "qseqid" |sort -k 1,1 -u > .table
awk 'BEGIN{OFS=FS="\t"} FNR==NR{a[$1]=$2 FS $3;next}{ print $0, a[$1]}' .table $ids |\
    sed 's/\t$//g' |\
    awk 'FS=OFS="\t" {if ( NF==2  ) print $0,"Missing"; else print $0}' > .newtable 
grep -Fwf $aed .newtable | cut -f 3 | sort | uniq -c | awk '{print $2"\t"$1}' > .counts
cut -f 3 .newtable |sort |uniq -c | awk '{print $2"\t"$1}' > .counts2
grep -Fwf $aed <(grep -vw "Zea_mays" .newtable | grep -wv "Tripsacinae" | grep -wv "Missing")  | cut -f 1,2 |sort | uniq > ids-pass.txt
grep -vw "Zea_mays" .newtable | grep -wv "Tripsacinae"| grep -wv "Missing" | cut -f 1,2 | sort | uniq > ids-And-pass.txt
awk 'BEGIN{OFS=FS="\t"} FNR==NR{a[$1]=$2 FS $3;next}{ print $0, a[$1]}' .order .counts2 | sort -k3,3 -n |cut -f 1,2 > $out2
awk 'BEGIN{OFS=FS="\t"} FNR==NR{a[$1]=$2 FS $3;next}{ print $0, a[$1]}' .order .counts | sort -k3,3 -n |cut -f 1,2 > $out

#rm .counts .ids .order .table .newtable
