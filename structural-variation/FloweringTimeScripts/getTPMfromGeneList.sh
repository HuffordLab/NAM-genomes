#!/bin/bash
head -n 1 B73-pasa-counts_tpm.txt | cut -f 1,6,7 > FTgenes_TPMcounts.txt
for i in *_tpm.txt; do
	grep -F -E -f grep-able_ids.txt ${i} | cut -f 1,6,7 > temp.txt ;
		awk -v OFS='\t' -v nam=${i%-pasa-counts_tpm.txt} '{print $0,nam}' temp.txt >> FTgenes_TPMcounts.txt ;
		echo Done with ${i} ;
	done 
		#grep -c -F -E ${line} ${i} ;
		#echo "transcript:Zm[0-9][0-9][0-9][0-9][0-9][a-z]${line:8}_T*" >> grep-able_ids.txt ;

TPM=$1
genome=${TPM%-pasa-counts_tpm.txt}
head -n 1 ${TPM} | cut -f 1,2,3,4,6,7 > ${genome}_candTPM.txt
slice=$(head -1 allgeneids_wB73v5.txt | tr ' ' '\n' | cat -n | grep ${genome} | cut -f 1)

tail +2 allgeneids_wB73v5.txt | cut -f ${slice} -d " " > temp_geneid.txt
awk 'BEGIN { FS=":";OFS="\n"} {print $1,$2,$3}' temp_geneid.txt |sed 's/NA//g' > temp_geneid2.txt
awk 'NF > 0 {print $0}' temp_geneid2.txt >temp_geneid3.txt
grep -F -E -f temp_geneid3.txt ${TPM} | cut -f 1,2,3,4,6,7 >> ${genome}_candTPM.txt
rm temp*.txt
echo Done with ${TPM}

cut -f 1 B97_candTPM.txt | cut -c 12-25 | head
grep Zm00018a002232 allgeneids_wB73v5.txt | cut -f 1 -d " "

TPM=$1
cut -f 1 ${TPM} | cut -c 12-25 > temp_gene.txt
grep -f temp_gene.txt allgeneids_wB73v5.txt | cut -f 1 -d " " >> temp_v3.txt
paste ${TPM} temp_v3.txt