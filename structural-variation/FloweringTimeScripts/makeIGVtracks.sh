#!/bin/bash
#make the IGV tracks from the SV bed files

#tail -n 1 V5_B73vsTzi8.bed | cut -f 15 | cut -f 11 -d : | awk -F - -v OFS=_ '{print $1,$2}'| awk -F _ '{print $1,$2,$4}'
#Gives the fields 1,2,3,7,8 of output track bed file
#import line of file | cut 15th field where the SV info is | cut to the 11th field, delimited by : |
#reprint that 11th field with _ as the delimiter | print out the SV chr start and stop info

#Field 4 (INS or DEL) get from 
#tail -n 1 V5_B73vsTzi8.bed | cut -f 15 | cut -f 7 -d : 

#Field 5 = .
#Field 6 get from tail -n 1 V5_B73vsTzi8.bed | cut -f 15 | cut -f 5 -d : 

#Field 9 if INS = 0,255,0 if DEL = 255,0,0
#while IFS=$'\t' read -r line
#do
#printf "%b\t" ${line} | cut -f 10 | cut -f 11 -d : | awk -F - -v OFS=_ '{print $1,$2}'| awk -F _ -v OFS='\t' '{print $1,$2,$4}' >> tmp_coords.bed
#printf "%b\t" ${line} | cut -f 10 | cut -f 7 -d : >> tmp_type.bed
#printf "%b\t" ${line} | cut -f 10 | cut -f 5 -d : >> tmp_strand.bed
#done < ${input}

#paste tmp_coords.bed tmp_type.bed tmp_strand.bed > tmp_track_${input}
#awk -v OFS='\t' '{print $1,$2,$3,$4,".",$5,$2,$3}' tmp_track_${input} > tmp_base_track_${input}
#echo made the base IGV track files

#awk -v OFS='\t' '$4 ~ /INS/ {print $0,"0,0,255"}' tmp_base_track_${input} >> tmp_color_track_${input} 
#awk -v OFS='\t' '$4 ~ /DEL/ {print $0,"255,0,0"}' tmp_base_track_${input} >> tmp_color_track_${input} 

#awk -v OFS='\t' '{gsub(/\,[a-zA-Z0-9_][a-zA-Z0-9_][a-zA-Z0-9_][a-zA-Z0-9_]/,"",$3)}1' tmp_color_track_${input} > tmp3.bed
#awk -v OFS='\t' '{gsub(/\,[a-zA-Z0-9_][a-zA-Z0-9_][a-zA-Z0-9_]/,"",$4)}1' tmp3.bed > tmp34.bed
#awk -v OFS='\t' '{gsub(/\,[a-zA-Z0-9_][a-zA-Z0-9_][a-zA-Z0-9_][a-zA-Z0-9_]/,"",$8)}1' tmp34.bed > color_track_${input} 

#echo finished making the color files with ${input}

#rm tmp*
#echo finished with ${input}, removed all tmp files 

#gen1=${input#*B73vs}
#gen=${gen1%.bed}

#echo track name="${gen}" description="Insertions or deletions in ${gen} relative to B73v5" > uniq_color_track_${input}
#sort color_track_${input} -k1,3 | uniq -u >> uniq_color_track_${input}
#echo Done with uniq_color_track_${input}

input=$1
module load bedtools2

#awk -v OFS='\t' '{print $2,$3,$4,$1}' ${input} > ${input%.bed}_wIDs.bed  #Adds a uniq ID to each row, orders for bedfile
#for i in *wIDs.bed ; do sort -k1,1 -k2,2n ${i} > sorted_${i} ; done #have to sort before bedtools merge
for i in *wIDs.bed; do awk -v OFS='\t' '{print $2,$3,$4,$1}' ${i} | sort -k1,1 -k2,2n > sorted_${i}; done #rearranges columns and sorts them
for i in sorted*.bed ; do bedtools merge -c 4 -o distinct -i ${i} > merged_${i} ; done #merging and printing the IDs
#cat extra_genes.bed promo_regions.bed gene_regions.bed | cut -f 1-4 > allregions.bed #to get all the regions of interest into 1 bed file
for i in merged_sorted_* ; do bedtools intersect -a allregions.bed -b ${i} -wa -wb > intersect_allregions_${i} ; done #doing the intersection

for i in intersect_allregions_merged_sorted_*INS*.bed ; do awk -v OFS='\t' '{print $5,$6,$7,$8,".",".",$6,$7,"0,0,255"}' ${i} > track_${i} ; done` #to make track files for intersected regions
for i in intersect_allregions_merged_sorted_*DEL*.bed ; do awk -v OFS='\t' '{print $5,$6,$7,$8,".",".",$6,$7,"255,0,0"}' ${i} > track_${i} ; done`
#for i in merged_sorted_*DEL*.bed; do awk -v OFS='\t' '{print $1,$2,$3,$4,".",".",$2,$3,"255,0,0"}' ${i} > track_${i} ; done` #make track files for the whole genome
#for i in merged_sorted_*INS*.bed; do awk -v OFS='\t' '{print $1,$2,$3,$4,".",".",$2,$3,"0,0,255"}' ${i} > track_${i} ; done`
#for i in B97 CML103 CML228 CML247 CML277 CML322 CML333 CML52 CML69 HP301 IL14H Ki11 Ki3 Ky21 M162W M37W Mo18W MS71 NC350 NC358 Oh43 Oh7b P39 Tx303 Tzi8 ; do cat track_merged_sorted_${i}*.bed > track_${i}.bed ; done` #putting Ins and Dels together
for i in B97 CML103 CML228 CML247 CML277 CML322 CML333 CML52 CML69 HP301 IL14H Ki11 Ki3 Ky21 M162W M37W Mo18W MS71 NC350 NC358 Oh43 Oh7b P39 Tx303 Tzi8 ; do for j in DEL DUP INS INV; do cat ${j}/track_intersect_allregions_merged_sorted_${j}_${i}_combined_wIDs.bed >> track_intersect_allregions_${i}.bed ; done; echo done with ${i};done
