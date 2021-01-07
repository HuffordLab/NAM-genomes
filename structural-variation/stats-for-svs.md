# file preprocessing 

``sed -e 's+AA+0/0+g' -e 's+TT+1/1+g' -e 's+NN+./.+g' NAM_founders_SVs.sniffles-bionano.hmp.txt > Nam_replaced.txt``
#### replace AA with 0/0(absent) and TT with 1/1(present) and NN with ./.(missing)

``cut -f 1 Nam_replaced.txt > SV_id.txt``
#### get SV_IDs

``cut -f 12- Nam_replaced.txt > Nam_genotypes.txt``
#### get SV_genotypes

``paste SV_id.txt Nam_genotypes.txt > combined.txt``
#### combine two files

``bash mkGenomeSV.sh combined.txt``
#### makeGenome SV files for each genome

``bash mkGenomeSV.sh combined.txt``
#### makeGenome SV files for each genome

``for i in *_combined.txt; do cut -f 1 ${i} > tmp_${i} ; done``
#### get rid of second column

``for i in tmp*; do tail -n +2 ${i} > tmp2_${i} ; done``
#### get rid of header

``for i in tmp2*; do sed 's/\./ /g' ${i} > ${i#tmp2_tmp_}.bed; done``
#### replace . with tab for SV_IDs

``for i in *txt.bed ; do mv ${i} ${i%.txt.bed}.bed ;done ``
#### changes the ending of the file names to .bed

``rm tmp*``
``rm  *_combined.txt``

####remove the temporary files

``for i in *.bed ; do awk -v OFS='\t' '{print $2,$3,$4,$1}' ${i} > ordered_${i} ; done``
#### reorder the bedfiles

``for i in ordered*; do grep tra ${i} >  TRA/tra_${i}; done``
``for i in ordered*; do grep -v tra ${i} > ORDERED_${i#ordered}; done``
``rm ordered*``
``for i in ORDERED*; do grep -v scaf ${i}> ordered_${i#ORDERED__}; done``
#### seperate out traslocations

``for i in ordered*; do sort -dsk1,1 -k2n,2 ${i} > sorted_${i} ; done``

#### sort the bedfiles before merge

``for i in sorted*; do bedtools merge -c 4 -o distinct -i ${i} > merged_${i} ; done``
#### merge bedfiles

``rm ordered*``
#### remove ordered_bedfiles

``rm sorted*``
#### remove sorted_bedfiles

``mkdir merged-files``
``mv merged* ./merged-files``
#### clean up directory

``mkdir INS`` ``mkdir DEL`` ``mkdir DUP`` ``mkdir INV``

``for i in merged*; do grep ins ${i} > ./INS/INS_${i} ; done``
``for i in INS*; do grep -v ',' ${i} > ins_${i#INS_}; done``
``rm INS*``

``for i in merged*; do grep del ${i} > ./DEL/DEL_${i} ; done``
``for i in DEL*; do grep -v ',' ${i} > del_${i#DEL_}; done``
``rm DEL*``

``for i in merged*; do grep dup ${i} > ./DUP/DUP_${i} ; done``
``for i in DUP*; do grep -v ',' ${i} > dup_${i#DUP_}; done``
``rm DUP*``

``for i in merged*; do grep inv ${i} > ./INV/INV_${i} ; done``
``for i in INV*; do grep -v ',' ${i} > inv_${i#INV_}; done``
``rm INV*``

## sv_count

 ``wc -l *bed > sv_count.bed``
 ``sed -e 's/ins_merged_sorted_ordered_//g' -e 's/_combined.bed//g' sv_count.bed >ins_count.bed``
 ``sed -e 's/inv_merged_sorted_ordered_//g' -e 's/_combined.bed//g' sv_count.bed >inv_count.bed``
 ``sed -e 's/del_merged_sorted_ordered_//g' -e 's/_combined.bed//g' sv_count.bed >del_count.bed``
 ``sed -e 's/dup_merged_sorted_ordered_//g' -e 's/_combined.bed//g' sv_count.bed >dup_count.bed``
#### counting numbers of svs for each line for each sv type
 ``for i in tra*; do awk '{print $1"\t"$2-5000"\t"$2+5000}' ${i} > tmp_${i}; done``
 ``for i in tmp*; do grep -v scaf ${i} > tmp2_${i}; done``
 ``for i in tmp2*; do sort -dsk1,1 -k2n,2 ${i} > sorted_${i#tmp2_tmp_} ; done``
 ``for i in sorted*; do bedtools merge -d 10000 -i ${i}> merged_${i}; done``
 ``rm tmp*`` ``rm tmp2*`` ``rm sorted*``
 `` wc -l *bed > sv_count.bed``
 ``sed -e 's/merged_sorted_tra_ordered_//g' -e 's/_combined.bed//g' sv_count.bed >tra_count.bed``
 *perform sv-t-test between tropical and temperate lines by sv-t-test.Rmd*
 
## sv_size

``for i in *combined.bed; do awk '{$5 = $3-$2} 1' ${i} > diff_${i} ;done``
#### Substract the stop and start position of the SV

``for i in diff*; do awk '{sum+=$5} END {print sum}' ${i} > sum_${i}; done``
#### calculating the sum of total sv_length for each sv type

``awk -v OFS='\t' '{print $0,FILENAME}' sum* > del_sum.bed``
``awk -v OFS='\t' '{print $0,FILENAME}' sum* > ins_sum.bed``
``awk -v OFS='\t' '{print $0,FILENAME}' sum* > inv_sum.bed``
``awk -v OFS='\t' '{print $0,FILENAME}' sum* > dup_sum.bed``

#### seperate different SVs

## Count by size
### 5kb and 500bp
``for i in merged*; do grep -v ',' ${i} > Merged_${i}; done``
``for i in Merged*; do awk '{$5 = $3-$2} 1' ${i} > diff_${i#Merged_} ;done``
``for i in diff*; do awk '$5<500' ${i} > 500bp/500bp_${i} ;done``
``for i in diff*; do awk '$5<5000' ${i} > 5kbp/5kb_${i} ;done``

## sv-gene-analysis

 *gene density in sv region*
1. get bed files(gene-model and sv) B73-pangene file 
2. use bedtools intersect
3. sum number of gene base pairs
4. sum number of sv basepairs
5. calculate genic basepair density by #gene bp sum/#sv bp sum
6. compare indel/inv/overall gene density


### Formatting files

``awk -v OFS='\t' '{print $2,$3,$4,$1}' curated_SVs_containing_genes.txt > svTMP``
``sed -i 's/^[^_]*_//' svTMP > sv.bed``

``awk '$5=$3-$2' sv.bed > diff_sv.bed``
``awk '$5>=1000000' diff_sv.bed > 1mb-sv.bed``
``cat less-than-1mb-sv.bed 1mb-sv.bed > full-sv.bed``
``grep -E 'DEL|del' full-sv.bed | sed 's/ /\t/g' | cut -f 1,2,3 | awk '$4=$3-$2' > DEL.bed``
``grep -E 'INV|inv' full-sv.bed | sed 's/ /\t/g' | cut -f 1,2,3 | awk '$4=$3-$2' > INV.bed``

``sed 's/ /\t/g' DEL.bed > del-final.bed``
``sed 's/ /\t/g' INV.bed > inv-final.bed``


### Intersect sv file and B73-pangene file

``bedtools intersect -wb -a del-final.bed -b B73-pangene.bed | sort | uniq > sv_genic_indel.bed``
``bedtools intersect -wb -a inv-final.bed -b B73-pangene.bed | sort | uniq > sv_genic_inv.bed``

``awk '$8=$3-$2' sv_genic_indel.bed | awk '{sum+=$8} END {print sum}' > gene-bp-del``
``awk '$8=$3-$2' sv_genic_inv.bed | awk '{sum+=$8} END {print sum}' > gene-bp-inv``


``awk '{sum+=$4} END {print sum}' del-final.bed > del-sum``
``awk '{sum+=$4} END {print sum}' inv-final.bed > inv-sum``

``awk '$4=$3-$2' B73-pangene.bed | awk '{sum+=$4} END {print sum}' > overall-gene-bp``


## Get number of large distinct svs(>1mb)
``awk '$5>=1000000' curated_SVs_containing_genes.txt > curated-1mb-sv.txt``
``awk -F"\t" '!seen[$3, $4]++' curated-1mb-sv.txt | grep INV > curated-1mb-inv.txt``
``awk -F"\t" '!seen[$3, $4]++' curated-1mb-sv.txt | grep DEL > curated-1mb-del.txt``


## MAF

### Variable denominator


``awk '$2+$3>0' temp > test ``
``awk -F '\t' '{if($2/($2+$3) < $3/($2+$3)) print $2/($2+$3); else print $3/($2+$3)}' test > tmp2``
``tail -n +2 test > filtered``
``awk '{split($1,a,";"); print a[3]}' filtered | sed 's/SVLEN=//g' | sed 's/-//g' > length``
``paste length tmp2 > maf-variable.txt``
#### reorder+filter

``awk '$1>100 && $1<=500' maf-variable.txt > 100bp-500bp``
``awk '$1>500 && $1<=5000' maf-variable.txt > 500bp-5000bp``
``awk '$1>5000 && $1<=10000' maf-variable.txt > 5kb-10kb``
``awk '$1>10000' maf-variable.txt > 10kb+``
#### bin maf values based on sv sizes

``awk -v OFS='\t' '{print $0,FILENAME}' 100bp-500bp > maf_binned.txt``
``awk -v OFS='\t' '{print $0,FILENAME}' 500bp-5000bp >> maf_binned.txt``
``awk -v OFS='\t' '{print $0,FILENAME}' 5kb-10kb >> maf_binned.txt``
``awk -v OFS='\t' '{print $0,FILENAME}' 10kb+ >> maf_binned.txt``
#### add size info

``cut -f 2,3 maf_binned.txt > MAF_binned.txt ``
``awk -v OFS='\t' '{print $2,$1}' MAF_binned.txt > maf_binned_v8.txt``
``awk '$2 > 0' maf_binned_v8.txt > maf_v8_variable.txt``
#### reorder
*generate plot by sv-maf-plot.Rmd*


## SV density


``genome.txt (for B73_V5)``

The genome file is based on B73_V5 for chromosome sizes


#### loading bedtools: 

``module load bedtools2``

#### make genome.windows bedfile from genome.txt file

``bedtools makewindows -b genome.txt -w 2500000 > genome.windows.bed``

for each 2500kb

-w is adjustable

#### make density bedfiles for each genome(cumulative)

``for i in merged*; do bedtools intersect -a genome.windows.bed -b ${i} -c > density_${i#merged_sorted_ordered_}; done``


#### make density bedfiles for each chromosome

``rename "density_" "" *``

``grep -w chr1 B97 CML247 CML333	HP301 Ki3 M37W NC350 Oh7b Tzi8 CML103 CML277 CML52 IL14H Ky21 MS71 NC358 P39 CML228 CML322 CML69 Ki11 M162W Mo18W Oh43 Tx303 > chr1.bed``
repeat for all chr1-chr10

``sed 's/:chr10//g' chr10.bed > final_10.bed``
repeat for all chr1.bed-chr10.bed

``for i in final*; do sed -e 's/Oh7b/Oh7B/g' -e 's/MS71/Ms71/g' -e 's/IL14H/Il14H/g' ${i} > ${i#final_};done``

*Generate plots by sv-density-plot.Rmd*
