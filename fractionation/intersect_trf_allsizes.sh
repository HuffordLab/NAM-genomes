#!/bin/bash

#This script compares the coordinates of the tandem array regions determined by trf with the coordinates of the filtered, formatted dc-megablast outputs and reports only non-intersecting dc-megablast alignment coordinates using bedtools intersect. The formatted dc-megablast output is first reformatted to a bed file, then formatted back to its original state after the bedtools step. This step could be optimized by formatting the dc-megablast output to a bedfile from the outset. 

module load bedtools2

for sample in *_ISUmasked_Sb_subgenomes_exons_dc-megablast_dagchainer.txt
        do
                echo $sample
                describer=$(echo ${sample} | sed 's/_ISUmasked_Sb_subgenomes_exons_dc-megablast_dagchainer.txt//')
                echo $describer


cat ${sample} | awk -v OFS="\t" '{print $5,$7,$8,$1,$2,$3,$4,$5,$9}' | bedtools sort | bedtools intersect -a - -b ${describer}*_trf_allsizes.bed -v | awk -v OFS="\t" '{print $4,$5,$6,$7,$1,$1"_"$2"_"$3,$2,$3,$9}' > ${describer}_ISUmasked_Sb_subgenomes_exons_dc-megablast_dagchainer_trf_allsizes_sb-trf.txt

done
