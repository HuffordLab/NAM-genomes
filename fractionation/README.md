# FRACTIONATION PIPELINE

**Key steps:**

1) Find tandem repeats in each masked genome using Tandem Repeat Finder
2) Align Sorghum exons associated with B73 subgenomes to each masked genome using dc-megablast
3) Remove dc-megablast alignments that overlap tandem repeat genomic coordinates as determined by Tandem Repeat Finder
4) Use the tandem-repeat filtered dc-megablast alignments from step 3) as inputs to DagChainer and run DagChainer
5) Organize and filter DagChainer outputs and generate tables


**Input files:**

a) Sbicolor_313_v3.1_exons_primary_notandems_cshl_clusters4.fa <br/>
b) Sb_exons_coords_CSHL_subgenomes_sections_sorted.txt<br/>
c) Repeatmasked NAM and B73 genomes<br/>


These steps were all done recursively on each masked NAM genome and B73.


**Software used for these analyses:**

a) Tandem Repeat Finder (TRF) 4.09 https://tandem.bu.edu/trf/trf.html<br/>
b) Blast+ 2.7.1<br/>
c) DagChainer http://dagchainer.sourceforge.net/<br/>
d) bedtools2 2.27 https://github.com/arq5x/bedtools2<br/>
e) R 3.6.0 (for post-analysis processing)<br/>
	- dplyr
	- purrr
	- tibble
	- tidyr
	- readr


## 1) Find Tandem Repeats


### Step 1: Run Tandem Repeat Finder (TRF)

```bash
sh trf.sh
```


**trf.sh** code

```bash
#!/bin/bash

for sample in *.fasta.masked
        do
                echo $sample
                describer=$(echo ${sample} | sed 's/.fasta.masked//')
                echo $describer

trf ${sample} 2 7 7 80 10 50 2000 -l 1 -d -h

done
```

### Step 2: Parse TRF output:


```bash
sh parse_trf_bedfile_allsizes.sh
```	

**parse_trf_bedfile_allsizes.sh** code

```bash
#!/bin/bash

for sample in *.fasta.masked.2.7.7.80.10.50.2000.dat
        do
                echo $sample
                describer=$(echo ${sample} | sed 's/.fasta.masked.2.7.7.80.10.50.2000.dat//')
                echo $describer

#Convert trf .dat file to .bed file using trf's TRFdat_to_bed.py script:

python TRFdat_to_bed.py --dat ${sample} --bed ${describer}.bed

#Parse above bed file:

cut -f 1,2,3 ${describer}.bed > ${describer}_trf_allsizes.bed

done
```




## 2) Align Sorghum exons to each masked genome using dc-megablast

### Step 1: Create blast database:

```bash
sh blastdb.sh
```

**blastdb.sh** code

```bash
#!/bin/bash

module load blast-plus

#looping blastdb:

for sample in *.fasta.masked
        do
                echo $sample
                describer=$(echo ${sample} | sed 's/.fasta.masked//')
                echo $describer

makeblastdb -in ${sample} -dbtype nucl -out ${describer}_masked

done
```


### Step 2: align reads using dc-megablast:


```bash
sh dc-megablast.sh
```

**dc-megablast.sh** code

```bash
#!/bin/bash

module load blast-plus

# have all your masked blast databases within the same directory

for database in *.nhr

        do

blastn -task dc-megablast -outfmt 6 -query Sbicolor_313_v3.1_exons_primary_notandems_cshl_clusters4.fa -db ${database%.*} -num_threads 4 -out "${database}_Sb3_exons_primary_notandems4_dc-megablast_nomax_ISUmasked.txt"

        done
```


### Step 3: parse dc-megablast output:

```bash
sh parse_dc-megablast.sh
```

**parse_dc-megablast.sh** code

```bash
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
```




## 3) Remove dc-megablast alignments that overlap Tandem Repeat Finder coordinates

```bash
sh intersect_trf_allsizes.sh
```


**intersect_trf_allsizes.sh** code

```bash
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
```



## 4) Use the tandem-repeat filtered dc-megablast alignments from step 3) as inputs to DagChainer and run DagChainer


```bash
sh dagchainer-filtered-trf.sh
```

**dagchainer-filtered-trf.sh** code

```bash
#!/bin/bash

for sample in *_ISUmasked_Sb_subgenomes_exons_dc-megablast_dagchainer_trf_allsizes_sb-trf.txt
        do
                echo $sample
                describer=$(echo ${sample} | sed 's/_ISUmasked_Sb_subgenomes_exons_dc-megablast_dagchainer_trf_allsizes_sb-trf.txt//')
                echo $describer


#First, do another filter for tandem repeats using Dagchainer's filter_repetitive_matches.pl script:

perl /DAGCHAINER/accessory_scripts/filter_repetitive_matches.pl 100000 < ${sample} > ${describer}_ISUmasked_Sb_subgenomes_exons_dc-megablast_dagchainer_trf_allsizes_sb-trf.filtered;

#Run dagchainer (the output will append the extension .aligncoords to each output file); parameters were optimized for maize, and -A 15 takes into consideration the fact that individual exons, not genes, are being processed:

perl /DAGCHAINER/run_DAG_chainer.pl -i ${describer}_ISUmasked_Sb_subgenomes_exons_dc-megablast_dagchainer_trf_allsizes_sb-trf.filtered -s -I -D 1000000 -g 40000 -A 15

        done
```



## 5) Organize and filter DagChainer outputs and generate table

This generates the exon count table of syntenic Sorghum aligned exons vs B73 and the NAMs within each syntenic block designated by DagChainer. Files were not sorted prior to processing since the syntenic block order from DagChainer was used. Exon counts were determined using bedtools groupby by grouping exons within a syntenic region that had the same Sorghum gene ID in that region and reporting the count of those exons, as well as the chromosome ID and the minimum coordinate (i.e. what results in the start coordinate) of those exons within the syntenic region. The result of this parsing was made into a table using R.


### Step 1: Get exon counts for each NAM genome and B73:

```bash
sh post_dagchainer_coords_csv_filtered_subg_unique_trf.sh
```

**post_dagchainer_coords_csv_filtered_subg_unique_trf.sh** code

```bash
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
```

### Step 2: Make table for each subgenome:

all M1 csv files were put into their own directory; same with all M2 csv files.

For each directory:

Module load r

```bash
./r_table.R
```

**r_table.R** code

```bash
#!/usr/bin/env Rscript

#make a table from all csv files in a directory

library(dplyr)
library(purrr)
library(tibble)
library(tidyr)
library(readr)

multmerge = function(mypath){filenames = list.files(path=mypath, full.names=TRUE)
datalist = lapply(filenames, function(x){read.csv(file=x,header=T)})
Reduce(function(x,y) {merge(x,y,all = TRUE)}, datalist)
}

full_data = multmerge("/path/to/directory")

write.csv(full_data, file = "output.csv")
```


### Step 3: Join subgenome tables:

Subgenome tables were joined together using unix join.
