#!/bin/bash
file=$1
ml bedtools2
snp1=/ptmp/LAS/arnstrm/sam/Permutation_Files/T2_sigSNPs.bed
snp2=/ptmp/LAS/arnstrm/sam/Permutation_Files/T30_sigSNPs.bed
snp3=/ptmp/LAS/arnstrm/sam/Permutation_Files/T3_sigSNPs.bed
bedtools intersect -a $file -b $snp1 -wa -wb  | cut -f 4 | sort | uniq -c | awk '{print $2"\t"$1}' > ${file%.*}_SNPs-intersecting-T2-gene-counts.txt
bedtools intersect -a $file -b $snp2 -wa -wb  | cut -f 4 | sort | uniq -c | awk '{print $2"\t"$1}' > ${file%.*}_SNPs-intersecting-T30-gene-counts.txt
bedtools intersect -a $file -b $snp3 -wa -wb  | cut -f 4 | sort | uniq -c | awk '{print $2"\t"$1}' > ${file%.*}_SNPs-intersecting-T3-gene-counts.txt
for f in T2 T30 T3; do
median=$(cut -f 2 ${file%.*}_SNPs-intersecting-${f}-gene-counts.txt | awk -f ~/bin/median.awk)
total=$(cut -f 4 $file | sort |uniq |wc -l)
hits=$(cut -f 1 ${file%.*}_SNPs-intersecting-${f}-gene-counts.txt |sort |uniq |wc -l)
echo -e "${file%.*}\t$f\t${total}\t${hits}\t${median}"
done
