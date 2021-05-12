#!/bin/bash
#goal is to make different bed files for each genome vs. B73
#need fields 1-14 and whichever field for each genome

#grep ^# NAM-structural-variations-v2.0.vcf | tail -n 1 #gets the header info off by 5 fields due to coordinate beds

counter=16 #has to be larger for intersect files because of the added columns for queryall
Del=$1 #deletion file
Ins=$2 #insertion file
vcf=$3 #vcf file with all the SV
nam=${Del%_*.bed}
ID=11 #starts at the column for B97
while [ ${counter} -le 40 ] #while the counter is less than 35 (end of the columns)
do
	gen=$(grep ^# ${vcf} | tail -n 1 | cut -f ${ID}) #cuts out only the genome of interest (a column)
	awk -v OFS='\t' -v c=${counter} '$c ~ /DEL/ {print}' ${Del} | 
	cut -f 1-9,${counter} >> ${nam}_B73vs${gen}.bed
	#awk finds all the deletions for that genome and then adds all the relevant info into a new bed file for that genome
counter=$(( ${counter} + 1 ))
ID=$(( ${ID} + 1 ))
echo Found all the DELs with ${gen}
done

counter=16
ID=11
while [ ${counter} -le 40 ]
do
	gen=$(grep ^# ${vcf} | tail -n 1 | cut -f ${ID})
	awk -v OFS='\t' -v c=${counter} '$c ~ /INS/ {print}' ${Ins} | 
	cut -f 1-9,${counter} >> ${nam}_B73vs${gen}.bed
counter=$(( ${counter} + 1 ))
ID=$(( ${ID} + 1 ))
echo Found all the INSs with ${gen}
done