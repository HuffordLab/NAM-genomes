#!/bin/bash
cds=$1
ml bioawk
bioawk -c fastx '{print $name"\t"substr($seq,0,3)"\t"reverse(substr(reverse($seq),0,3))}' ${cds} > ${cds}-start-stop.txt
sed -i 's/_/\t/1' ${cds}-start-stop.txt
awk '$3  ~ /[ATGC]TG/ && $4  ~ /TAA|TAG|TGA/' ${cds}-start-stop.txt | sort -k1,1 -u | sed 's/\t/_/1' > ${cds}-reassigned-transcript.txt
gnum=$(cut -f 1 ${cds}-start-stop.txt | sort | uniq |wc -l);
onum=$(awk  '$3  ~ /[ATGC]TG/ && $4  ~ /TAA|TAG|TGA/'  ${cds}-start-stop.txt | cut -f 1| sort | uniq | wc -l);
nam=$(echo $cds | cut -f 2 -d "-")
echo -e "${nam}\t${gnum}\t${onum}" | awk '{print $0"\t"$3/$2}'
