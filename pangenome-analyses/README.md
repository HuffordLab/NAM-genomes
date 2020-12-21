# NAM_Pan_Genome_Analysis

This file includes the method used for maize pan genome construction using 26 NAM genome asseblies.

The method includes the following steps: 
```
1. using nucmer identify sytenic block
2. using split-merge pipeline identify one to one gene mapping and tandem duplicates 
3. using R join the gene linking information and merge/merge genes information and construct pan gene matrix 
```

Module involved in the pan-genome construction:
```
nucmer: mummer/4.0.0.beta2
split_merge pipeline: https://github.com/HirschLabUMN/Split_genes/tree/master/Split_Merge_Pipeline
ncbi_blast: ncbi_blast+/2.8.1
```

Note:
```
Genes used for this pan genome analysis are genes annoated on chromsome, scaffold genes are not included. 
```

Pan-genome construction:
1. Extract canonical sequence ID that are on chromsomes from the gff file 
```
for i in *.gff; do
  grep "canonical_transcript" "$i" | grep -v scaf |cut -f 9 | tr ";" "\n" | grep "transcript_id" | sed 's/transcript_id=//g' > "$i"_canonical_transcript.txt
done  
```
2. Pulling canonical transcript gene gff file, below is one example 
```
perl pull_out_caononical_transcript_coded_from_gff.pl -i ~/nam_pan_genome/NAM_annotation/gff/zea_maysb73_core_3_87_1.gff -l  ~/nam_pan_genome/NAM_annotation/canonical_transcript/B73_canonical_transcript.txt -o B73_canonical.gff

scaf transcript were removed via grep -v scaf * 
```
3. Using canonical transcript CDS sequences to make blast database 
```
perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' *_canonical_transcript_ID zea_mays_*_core_3_87_1.canonical.cds.fasta >  cds fasta without scaf

makeblastdb -in zea_maysb73_core_3_87_1.canonical.cds.fasta -out B73_cds_db -dbtype nucl
```
4. running nucmer at -c 1000, example of the script is shown below. All pairs are generated using mummer_c1000_propogator.sh 
```
nucmer --mum -c 1000 -p HP301_B73_c1000 /panfs/roc/groups/6/maize/shared/databases/genomes/Zea_mays/HP301_NAMassembly/Zm-HP301-REFERENCE-NAM-1.0/Zm-HP301-REFERENCE-NAM-1.0.fasta /panfs/roc/groups/6/maize/shared/databases/genomes/Zea_mays/B73_NAMassembly/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.fasta
```
nucmer output was further processed using nucmer_post_processing.sh
and filtered that the syntenic block must be on the same chr
```
for j in *.coords ; do
  grep -v scaf "$j"  | awk '(NR>1) && ($12 == $13 )' > "$j"_filter_chr_nucmer 
done 
```

5. Set up all by all blast to search pan gene using script all_by_all_blast_batch.sh. Below is an example of the script
```
python3 All_by_All_Blast_COedits_10.py -q /home/hirschc1/qiuxx221/nam_pan_genome/NAM_annotation/canonical_gff/B73_canonical.gff -s /home/hirschc1/qiuxx221/nam_pan_genome/NAM_annotation/canonical_gff/P39_canonical.gff -b /home/hirschc1/qiuxx221/nam_pan_genome/NAM_annotation/canonical_fasta/P39_cds_db -l /home/hirschc1/qiuxx221/nam_pan_genome/NAM_annotation/canonical_fasta/zea_maysb73_core_3_87_1.canonical.cds.fasta -o /home/hirschc1/qiuxx221/nam_pan_genome/all_by_all_1000_new_annotation_10_filter_nucmer/B73_P39_AllbyAll_res.txt -n /home/hirschc1/qiuxx221/nucmer_1000_filter/B73_P39_c1000.fil.coords -g ab
```
After running the first 4 modules, output is reshaped for R processing

6. Add header to each of the canonical transcript file so the list can be used as query list 
```
sed -i '1 i\Query_gene' *_canonical_transcript.txt 
```
7. Duplicate the canonical transcript ID so it can be used for left_join in R
```
for i in *.txt; do
  id=$(echo "$i" | cut -d'_' -f1)
  paste -d '\t' "$i" "$i"| sed "1iQuery_gene\t$id" >  "$i"_add.txt
done 
```

8. processing all_by_all_blast output 
```
# optional, could apply filter on the one to one mapping type
for j in *.txt ; do
  grep 'adjacent_genes_syntenic' "$j" | cut -f 1,12 | awk '{split($2,a,";");for(i in a)print $1"\t"a[i]}' | awk -F"," '$1=$1' OFS="\t" | cut -f 1,2 | sed 's/\t/,/g' >  $(basename "$j")_adj_pairs.txt 
  awk -F":" '$1=$1' OFS="\t" "$j" | cut -f 1-18 | grep -v NA | awk 'NR >1, $19=$18/$16'| awk '(NR>0) && ($19 >= 0.0 )'| sed 's/ /,/g' | grep one_to_one_mapping | cut -d ',' -f 1,12 >  $(basename "$j")_flt121_pair.txt 
  cat $(basename "$j")_adj_pairs.txt $(basename "$j")_flt121_pair.txt  > $(basename "$j")_combine.txt
done

# carry out the files named as _combine.txt, this file has pairs from both adjacent and one to one mapping type

# change comma separated rows into "\t" separated rows 
sed -i 's/,/\t/g' *

# sort all gene ID before the first merging by reference 

for j in *.txt; do
 sort "$j" > $(basename "$j")_sort.txt
done 

# for duplicates, merge back to semi comma separated 
for j in *sort.txt; do
  awk '
      BEGIN{FS="\t"; OFS=FS}; 
      { arr[$1] = arr[$1] == ""? $2 : arr[$1] ";" $2 }   
      END {for (i in arr) print i, arr[i] }
    ' "$j" > $(basename "$j")_merge.txt
done 

# add header to the gene pairs before downloading for R processing 
for i in *.txt; do
    id=$(echo "$i" | cut -d'_' -f2)
    sed "1iQuery_gene\t$id" $i > $(basename "$i")_fmt_col_name.txt
done

# rename file 
rename _AllbyAll_res.txt_combine.txt_sort.txt_merge.txt_fmt_col_name '_all_pair' *.txt
```

9. input all output from step 6,7,8 into R to make the base matrix using R script: make_base_pan_matrix.R
```
This step input all pairs information and align into a single matrix. gene without pairs is recored as NA 
```

10. compressing the base matrix by merging ID that present more than once within a column using R script: compressing_duplicate_ID.R
```
Gene pairs will be picked up by all by all blast files in both directions (B73_B97 vs B97_B73). gene ID has the chance to present more than once. This step joins all gene ID, tandem dulicates will be separated by ";"
```

11. for genes with NA value, further step is done to validate if it does exist in the genomic DNA (scaf_free) via gmap
```
#extract pan_gene_ID using matrix from step 10
cut -d "," -f 2 pan26_all.collapsed.csv | cut -d ";" -f 1 > pan_gene_ID.txt
#extract fasta sequence using pan_gene_id. The master_26.fasta is a concatenate file from catting all genomecds.fasta 
perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' pan_gene_ID.txt /home/hirschc1/qiuxx221/nam_pan_genome/NAM_annotation/canonical_fasta/master_26.fasta >  master_pan.fasta 

#using gmap makde database for each nam genome. cmd example as below
gmap_build -D /home/hirschc1/qiuxx221/nam_pan_genome/gmap_db -d B73 Zm-B73-REFERENCE-NAM-5.0_1.fasta

#map pan_fasta to all indivudual NAM genomes and output gff3 with a single path for each alginment. below is an example of the cmd 
gmap -D /home/hirschc1/qiuxx221/nam_pan_genome/gmap_db -d B73 -n 1 -f gff3_gene master_pan.fasta >  pan_to_B73_1_path.gff3
```

12.  output the alignment and apply filter on gmap coverage and identity was done using script output_gmap_canonical_gene_position.py 
```
#output format: 
pan_gene_ID NAM gmap_ID= genome coordinate 
#an gff format output was also generated to intersect with existing gene model using this script 
```

13. When the gmap CDS is larger than 200 bp and it can be intersected with NAM annotation gff CDS, the gmap coordinates were changed back to the gene name.

```
Refer to step two in gmap_pipeline section
```

14. the intersected gff files was further filtered using script intersect_gene_model_name_changing.py
```
This script changed coordinates back to NAM gene ID, if the mapping coverage covered at least 90% of the gene model

The final output format will be 
Pan_gene_ID NAM gene coordiates(either gmap_ID=xxx or Zm00xxxx)
```
15. convert output from step 14 using R script rejoin_matrix.R, but this is modified version that only allow to merge the pan gene ID. there is no compression for any possible duplicate gene IDs. They are removed manually when all the gene ID are identical
```
This step generate the final matrix for the pan genome that includes NA genes been recovered using gmap approach. A list of manually removed pan genes can be find in data visualization under duplicated_lines_manually_removed.csv file
```

16. Note on subgenome information
```
When one pan gene has more than 1 hits with different subgenome outcome, the majority subgenome genome is kept. For example, when pan gene A has 2 hits to M1, 1 hit to M2, M1 is considered as the subgenome where the pan gene is coming from. When M1, and M2 has the same hit, M1 is kept. Only 42 pan gene has this issue. The details of the conflicting subgenome information can be found at subgenome_conflict_select.csv
```
