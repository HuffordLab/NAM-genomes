#!/bin/bash
file=$1
#This turns txt file into bed file, columns = chromosome, start, end, genenames, strand
#echo ${file} | tr '\r' '\n' | #trims out headers and weird end characters
awk -v OFS="\t" '{print $3,$4,$5,$2,$6}'| #prints out fields in desired order
awk '{gsub(/chr|\"/,""); print}' | #removes chr from chromosome column and "" from rest of file
sort -n -k1,1 -k2,2 > ${file%.txt}.bed #sorts file and writes to bed file name

#This grabs the 5kb upstream of the genes from the bed file

head -n 1 ${file%.txt}.bed >> promoter_${file%.txt}.bed
grep + ${file%.txt}.bed > tmp.bed
awk -v OFS='\t' '{print $1, $2-5000, $2,$4,$5}' tmp.bed >> promoter_${file%.txt}.bed
echo Finished with positive strand
rm tmp.bed
grep - ${file%.txt}.bed > tmp.bed
awk -v OFS='\t' '{print $1, $2+5000, $2,$4,$5}' tmp.bed >> promoter_${file%.txt}.bed
echo Finished with negative strand
rm tmp.bed
	