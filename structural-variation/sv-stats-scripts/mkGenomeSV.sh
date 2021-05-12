#!/bin/bash

input=$1
counter=2
echo Starting with making individual genome SV files with ${input}
while [ ${counter} -le 27 ] 
	do
	gen=$(head -n 1 ${input} | cut -f ${counter})
	tail -n +2 ${input} | cut -f 1,${counter} > tmp_${gen}.txt
	head -n 1 ${input} | cut -f 1,${counter} > ${gen}_${input}
	grep "1\/1" tmp_${gen}.txt >> ${gen}_${input}
counter=$(( ${counter} + 1 ))
echo Done with ${input} for ${gen}
done
rm tmp*.txt
echo Done making individual genome SV files with ${input}
