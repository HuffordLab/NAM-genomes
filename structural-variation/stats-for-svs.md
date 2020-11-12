
# 1. File Processing 
*input File: NAM_founders_SVs.sniffles-bionano.hmp.txt*

``sed -e 's+AA+0/0+g' -e 's+TT+1/1+g' -e 's+NN+./.+g' NAM_founders_SVs.sniffles-bionano.hmp.txt > Nam_replaced.txt``
#### replace AA with 0/0(absent) and TT with 1/1(present) and NN with ./.(missing)

``cut -f 1 Nam_replaced.txt > SV_id.txt``
#### get SV_IDs

``cut -f 12- Nam_replaced.txt > Nam_genotypes.txt``
#### get SV_genotypes

``paste SV_id.txt Nam_genotypes.txt > combined.txt``
#### combine two files

``awk -v OFS='\t' '{print NR-1 $0}' combined.txt > combined2.txt ``
#### add id number to the beginning of the first field

``grep -v tra combined2.txt > combined_wIDs.txt``
#### exclude tra

``bash mkGenomeSV.sh combined_wIDs.txt``
#### makeGenome SV files for each genome

``for i in *_combined_wIDs.txt; do cut -f 1 ${i} > tmp_${i} ; done``
#### get rid of second column

``for i in tmp*; do tail -n +2 ${i} > tmp2_${i} ; done``
#### get rid of header

``for i in tmp2*; do sed 's/\./ /g' ${i} > ${i#tmp2_tmp_}.bed; done``
#### replace . with tab for SV_IDs

``for i in *txt.bed ; do mv ${i} ${i%.txt.bed}.bed ;done ``
#### changes the ending of the file names to .bed

``rm tmp*``
``rm *wIDs.txt``

#### remove the temporary files

``for i in *wIDs.bed ; do awk -v OFS='\t' '{print $2,$3,$4,$1}' ${i} > ordered_${i} ; done``
#### reorder the bedfiles

``for i in ordered*; do sort -dsk1,1 -k2n,2 ${i} > sorted_${i} ; done``
#### sort the bedfiles before merge

``for i in sorted_ordered*; do bedtools merge -c 4 -o distinct -i ${i} > merged_${i} ; done``
#### merge bedfiles

``rm ordered*``
#### remove ordered_bedfiles

``rm sorted*``
#### remove sorted_bedfiles

``mkdir merged-files``
``mv merged* ./merged-fles``
#### clean up directory

``mkdir INS`` ``mkdir DEL`` ``mkdir DUP`` ``mkdir INV``
``for i in *wIDs.bed; do grep ins ${i} > ./INS/INS_${i} ; done``
``for i in *wIDs.bed; do grep del ${i} > ./DEL/DEL_${i} ; done``
``for i in *wIDs.bed; do grep dup ${i} > ./DUP/DUP_${i} ; done``
``for i in *wIDs.bed; do grep inv ${i} > ./INV/INV_${i} ; done``
#### seperate different SVs


# 2. SV Count Analysis

 ``wc -l *bed > sv_count.bed``
 ``sed -e 's/INS_merged_sorted_ordered_//g' -e 's/_combined_wIDs.bed//g' sv_count.bed >ins_count.bed``
 ``sed -e 's/INV_merged_sorted_ordered_//g' -e 's/_combined_wIDs.bed//g' sv_count.bed >inv_count.bed``
 ``sed -e 's/DEL_merged_sorted_ordered_//g' -e 's/_combined_wIDs.bed//g' sv_count.bed >del_count.bed``
 ``sed -e 's/DUP_merged_sorted_ordered_//g' -e 's/_combined_wIDs.bed//g' sv_count.bed >dup_count.bed``
#### counting numbers of svs for each line for each sv type
run `sv_summary_stat_v8.Rmd` for corresponding plots



# 3. SV Size Analysis

``for i in *wIDs.bed; do awk '{$5 = $3-$2} 1' ${i} > diff_${i} ;done``
#### Substract the stop and start position of the SV

``for i in diff*; do awk '{sum+=$5} END {print sum}' ${i} > sum_${i}; done``
#### calculating the sum of total sv_length for each sv type

``awk -v OFS='\t' '{print $0,FILENAME}' sum* > del_sum.bed``
``awk -v OFS='\t' '{print $0,FILENAME}' sum* > ins_sum.bed``
``awk -v OFS='\t' '{print $0,FILENAME}' sum* > inv_sum.bed``
``awk -v OFS='\t' '{print $0,FILENAME}' sum* > dup_sum.bed``

run `sv_summary_stat_v8.Rmd` for corresponding plots


# 4. Cumulative SV count analysis(bootstrap500)

``grep -v tra combined.txt > combined2.txt``
#### exclude tra

``bash mkGenomeSV.sh combined_wIDs.txt``
#### makeGenome SV files for each genome

``for i in *_combined2.txt; do cut -f 1 ${i} > tmp_${i} ; done``
#### get rid of second column

``for i in tmp*; do tail -n +2 ${i} > tmp2_${i} ; done``
#### get rid of header


``for k in {1..50}; do for j in {1..10}; do for i in $(ls *.txt | shuf); do cat $i >> temp.$j.$k; sort -u temp.$j.$k|wc -l; done|perl transpose3.pl - > result.$j.$k; rm temp.$j.$k; done & done``

#### bootrap500times

``cat result* > sv_bootstrap500.txt``
run `bootstrapSV_v8.Rmd` for corresponding plots


# 5. MAF Analysis

``cut -f 2 site-summary.txt| sed 's/\./\t/g' > temp.txt``
``cut -f 3,4 temp.txt | tail -n +2 > temp1.txt``
``awk '{$3 = $2-$1} 1' temp1.txt | cat | cut -d " " -f 3 > temp2.txt``
``cut -f 15 site-summary.txt | tail -n +2 > temp3.txt``
``paste temp2.txt temp3.txt > maf_total``
#### formatting

``awk '$1>100 && $1<=500' maf_total > 100bp-500bp``
``awk '$1>500 && $1<=5000' maf_total > 500bp-5000bp``
``awk '$1>5000 && $1<=10000' maf_total > 5kb-10kb``
``awk '$1>10000' maf_total > 10kb+``
#### bin maf values based on sv sizes

``awk -v OFS='\t' '{print $0,FILENAME}' 100bp-500bp > maf_binned.txt``
``awk -v OFS='\t' '{print $0,FILENAME}' 500bp-5000bp >> maf_binned.txt``
``awk -v OFS='\t' '{print $0,FILENAME}' 5kb-10kb >> maf_binned.txt``
``awk -v OFS='\t' '{print $0,FILENAME}' 10kb+ >> maf_binned.txt``
#### add size info

``cut -f maf_binned.txt > MAF_binned.txt``
``awk -v OFS='\t' '{print $2,$1}' MAF_binned.txt > maf_binned_v8.txt``
``awk '$2 > 0' maf_binned_v8.txt > maf_v8.txt``
#### reorder
run `bootstrapSV_v8.Rmd` for corresponding plots


# 6. SV Density Analysis

``genome.txt (for B73_V5)``

The genome file is based on B73_V5 for chromosome sizes

#### loading bedtools: 

``module load bedtools2``

#### make genome.windows bedfile from genome.txt file

``bedtools makewindows -b genome.txt -w 2500000 > genome.windows.bed``

for each 2500kb

-w is adjustable

#### make density bedfiles for each genome(cumulative)

``for i in *bed; do bedtools intersect -a genome.windows.bed -b ${i} -c > density_${i#merged_sorted_ordered_}; done``

#### make density bedfiles for DEL of each genome

``for i in *bed; do grep del ${i} > DEL_${i#merged_sorted_ordered_}; done``
``for i in *bed; do bedtools intersect -a genome.windows.bed -b ${i} -c > density_${i}; done``

#### make density bedfiles for INS of each genome

``for i in *bed; do grep ins ${i} > INS_${i#merged_sorted_ordered_}; done``
``for i in *bed; do bedtools intersect -a genome.windows.bed -b ${i} -c > density_${i}; done``

#### make density bedfiles for each chromosome
``rename "density_" "" *``
``rename "_combined_wIDs.bed" "" *``

``grep -w chr1 B97 CML247 CML333	HP301 Ki3 M37W NC350 Oh7b Tzi8 CML103 CML277 CML52 IL14H Ky21 MS71 NC358 P39 CML228 CML322 CML69 Ki11 M162W Mo18W Oh43 Tx303 > chr1.bed``
repeat for all chr1-chr10

``sed 's/:chr1//g' chr1.bed > final_1.bed``
repeat for all chr1.bed-chr10.bed

run `sv_density_v8.Rmd` for corresponding plots
