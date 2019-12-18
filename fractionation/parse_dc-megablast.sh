#!/bin/bash

#This script formats the dc-megablast outfmt 6 output and merges it with the coordinates associated with the Sorghum exons based on Sorghum exon ID for preparation for DagChainer; DagChainer needs a very specific format, as specified in their documentation. Sorghum exons not on chromosomes are filtered out. The Sorghum exons in Sb_exons_coords_CSHL_subgenomes_sections_sorted.txt have the subgenome and chromosome information vs B73 associated with them; alignments to the NAMs to Sorghum that don't match the B73 chr ID were ignored. 

module load bedtools2

for sample in *_Sb3_exons_primary_notandems4_dc-megablast_nomax_ISUmasked.txt
        do
                echo $sample
                describer=$(echo ${sample} | sed 's/_Sb3_exons_primary_notandems4_dc-megablast_nomax_ISUmasked.txt//')
                echo $describer

cat ${sample} | grep -v "Sobic.K" | grep -v "super" | awk -v OFS="\t" '{if($9>$10) print $1,$2,$10,$9,$11; else print $1,$2,$9,$10,$11}' | sort -k1,1 | join - Sb_exons_coords_CSHL_subgenomes_sections_sorted.txt | tr ' ' '\t' | awk -v OFS="\t" '{if($2==$11) print $6,$1"_"$9"_"$10"_"$11,$7,$8,$2,$2"_"$3"_"$4,$3,$4,$5}' > ${describer}_ISUmasked_Sb_subgenomes_exons_dc-megablast_dagchainer.txt

	done
