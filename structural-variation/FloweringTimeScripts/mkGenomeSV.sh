#!/bin/bash

input=$1
counter=7
echo Starting with making individual genome SV files with ${input}
while [ ${counter} -le 32 ] 
	do
	gen=$(head -n 1 ${input} | cut -f ${counter})
	tail -n +2 ${input} | cut -f 1-6,${counter} > tmp_${gen}.txt
	head -n 1 ${input} | cut -f 1-6,${counter} > ${gen}_${input}
	grep "1\/1" tmp_${gen}.txt >> ${gen}_${input}
counter=$(( ${counter} + 1 ))
echo Done with ${input} for ${gen}
done
rm tmp*.txt
echo Done making individual genome SV files with ${input}

echo Starting making bed files
for k in *_${input}; do
	while IFS=$'\t' read -r line 
		do
		printf "%b\t" ${line} | cut -f 2 | cut -f 1 -d : | cut -f 1 -d - >> tmp_chr.txt
		printf "%b\t" ${line} | cut -f 2 | cut -f 1 -d : | cut -f 2 -d - >> tmp_start.txt
		printf "%b\t" ${line} | cut -f 2 | cut -f 2 -d : | cut -f 2 -d - >> tmp_end.txt
		printf "%b\t" ${line} | cut -f 1 >> tmp_id.txt
		done < ${k}
	paste tmp_chr.txt tmp_start.txt tmp_end.txt tmp_id.txt >> ${k%.txt}.bed
	rm tmp*.txt
	echo Done making bed files for ${k}
done
echo Done making all bed files