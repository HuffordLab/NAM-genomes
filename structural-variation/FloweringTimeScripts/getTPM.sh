#!/bin/bash

TPM=$1
genome=${TPM%-counts_tpm.txt}
sed 's/_R1_round.2Aligned.sortedByCoord.out.bam//g' ${TPM} | head -n 1 > ${genome}_candTPM.txt
slice=$(head -1 allgeneids_wB73v5.txt | tr ' ' '\n' | cat -n | grep ${genome} | cut -f 1)

tail -n +2 allgeneids_wB73v5.txt | cut -f ${slice} -d " " > temp_geneid.txt
awk 'BEGIN { FS=":";OFS="\n"} {print $1,$2,$3}' temp_geneid.txt |sed 's/NA//g' > temp_geneid2.txt
awk 'NF > 0 {print $0}' temp_geneid2.txt >temp_geneid3.txt
grep -F -E -f temp_geneid3.txt ${TPM} >> ${genome}_candTPM.txt
rm temp*.txt
echo Done with ${TPM}

#### With the read coverage files

TPM=$1
genome=${TPM%_joined_cnts_wLengths_tpm.txt}
head -n 1 ${TPM} > ${genome}_candTPM.txt

#This gets the column with the appropriate gene names 
slice=$(head -1 allgeneids_wB73v5.txt | tr ' ' '\n' | cat -n | grep ${genome} | cut -f 1)

tail -n +2 allgeneids_wB73v5.txt | cut -f ${slice} -d " " > temp_geneid.txt
awk 'BEGIN { FS=":";OFS="\n"} {print $1,$2,$3}' temp_geneid.txt |sed 's/NA//g' > temp_geneid2.txt
awk 'NF > 0 {print $0}' temp_geneid2.txt >temp_geneid3.txt
grep -F -E -f temp_geneid3.txt ${TPM} >> ${genome}_candTPM.txt
rm temp*.txt
echo Done with ${TPM}