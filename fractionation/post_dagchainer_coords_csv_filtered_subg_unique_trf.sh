#!/bin/bash

#Each subgenome was processed individually. Data was converted to csv files. 

module load bedtools2

for sample in *_ISUmasked_Sb_subgenomes_exons_dc-megablast_dagchainer_trf_allsizes_sb-trf.filtered.aligncoords
        do
                echo $sample
                describer=$(echo ${sample} | sed 's/_ISUmasked_Sb_subgenomes_exons_dc-megablast_dagchainer_trf_allsizes_sb-trf.filtered.aligncoords//')
                echo $describer

cat ${sample} | grep -v "#" | grep "M1" | awk '{$1=$1} 1' FS=".1.v3.1." OFS="\t" | tr '_' '\t' | awk -v OFS="\t" '{print $2".1_"$4"_"$5,$10,$11}' | groupBy -i - -g 1,2 -c 1,3 -o count,min  | sort -k1,1 -k3,3r | awk '!x[$1]++' | awk -v x=$describer 'BEGIN { OFS="\t"; print "gene_chr","pos_"x,"exoncount_"x} { print $1"_"$2,$4,$3}' | tr '\t' ',' > ${describer}_ISUmasked_Sb_dagchainer_filtered_dc-m_coords_trf_M1.csv

cat ${sample} | grep -v "#" | grep "M2" | awk '{$1=$1} 1' FS=".1.v3.1." OFS="\t" | tr '_' '\t' | awk -v OFS="\t" '{print $2".1_"$4"_"$5,$10,$11}' | groupBy -i - -g 1,2 -c 1,3 -o count,min  | sort -k1,1 -k3,3r | awk '!x[$1]++' | awk -v x=$describer 'BEGIN { OFS="\t"; print "gene_chr","pos_"x,"exoncount_"x} { print $1"_"$2,$4,$3}' | tr '\t' ',' > ${describer}_ISUmasked_Sb_dagchainer_filtered_dc-m_coords_trf_M2.csv

done
