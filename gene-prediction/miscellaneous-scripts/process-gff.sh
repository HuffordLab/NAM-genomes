#!/bin/bash
ml bedops
gff=$1
genome=$2
gff2bed < $gff > ${gff%.*}.bed
#awk -v x=${gff%.*} 'BEGIN{OFS=FS="\t"}{print>$8"."x".bed"}' ${gff%.*}.bed
awk 'BEGIN{OFS=FS="\t"} $8=="five_prime_UTR"' ${gff%.*}.bed | awk '{OFS=FS="\t"}$4=$10' > ${gff%.*}_5UTR.bed
awk 'BEGIN{OFS=FS="\t"} $8=="three_prime_UTR"' ${gff%.*}.bed | awk '{OFS=FS="\t"}$4=$10' > ${gff%.*}_3UTR.bed
ml purge
ml bedtools2
bedtools getfasta -fi $genome -fo ${gff%.*}_5UTR.fa -bed ${gff%.*}_5UTR.bed -name -s
bedtools getfasta -fi $genome -fo ${gff%.*}_3UTR.fa -bed ${gff%.*}_3UTR.bed -name -s
