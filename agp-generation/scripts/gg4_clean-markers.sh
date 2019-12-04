#!/bin/bash
echo "note that this will randomly select new set of markers everytime you run"
rm  pg-markers-mapped.csv &> /dev/null
marker=$1
nam=$(echo $marker |sed 's/_GG-mapped.csv//g')
grep -v "Scaffold" ${marker}  | cut -f 1 -d "," | sort | uniq -c | awk '{print $2"\t"$1}' | sort -k1,1 -n > ${nam}-marker-density.txt
awk 'BEGIN{OFS=FS="\t"}FNR==NR{a[$1]=$2 FS $3;next}{ print $0, a[$1]}' ../${nam}_star/${nam}.scaffolds.fasta.sizes ${nam}-marker-density.txt |awk '{print $0"\t"($2/$3)*10000}' > ${nam}-marker-table.txt
awk '$3>100000 && $2>20 && $3>1' ${nam}-marker-table.txt | cut -f 1 > scafs-for-agp.txt

while read line;  do
   awk -v x=$line -F"," '$1==x'  ${nam}_GG-mapped.csv|sort |uniq | sort -R | head -n 100;
done<scafs-for-agp.txt >> pangenome.csv
cp ../GoldenGate/*_GG-mapped.csv ./goldengate.csv
