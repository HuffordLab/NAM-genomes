#!/bin/bash
file=$1
#mark start as yes or no
awk 'FS=OFS="\t" {if ( $3~/.TG/ ) print $1,$2,"yes",$4; else print $1,$2,"no",$4}' $file > $file.temp1
awk 'FS=OFS="\t" {if ( $4~/TA[AG]/ || $4~/TGA/  ) print $1,$2,$3,"yes"; else print $1,$2,$3,"no"}' $file.temp1 > $file.temp2
awk 'FS=OFS="\t" {if ( $3=="yes" && $4="yes"  ) print $0,"yes"; else print $0,"no"}' $file.temp2 > ${file%.*}.binary.txt
rm $file.temp1 $file.temp2
exon=$(echo $file |sed 's/.len-start-stop.txt/-exon.bed.counts/g')
awk 'BEGIN{OFS=FS="\t"}FNR==NR{a[$1]=$2 FS $3;next}{ print $0, a[$1]}' $exon ${file%.*}.binary.txt | sed 's/\t$//g' > ${file%.*}.binary-with-exon-counts.txt

