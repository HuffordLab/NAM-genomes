---
title: "NAM conservation analysis readme"
author: "Asher Hudson"
date: "7/20/2020"
output:
  html_document: default
  pdf_document: default
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, eval = FALSE)
```

# Conservation analysis pipeline

### Key steps

1. Align query genomes to maize B73 v5 reference genome.
2. Use GERP to identify conserved elements.
3. Statistical analysis to detect enrichment/depletion of structural variants in conserved elements.


### Software used for these analyses *(add versions)*

a) Last (lastdb, last-train, lastal, maf-convert, last-split, maf-swap, last-postmask) v 2.31.1
b) Kent-utils (axtChain, chainMergeSort, faSize, chainPreNet, chainNet, faToTwoBit, netToAxt, axtToMaf)
c) multiz (multiz.v10.6)
d) GERP++
e) bedtools 2.27.1
f) R 3.6.3
g) R packages
    i) rPHAST_1.6.9
    ii) dplyr_0.8.5
    iii) data.table_1.12.6
    iv) RColorBrewer_1.1-2
    v) tidyr_1.0.0
    vi) ggplot2_3.3.0
    vii) grid_3.5.1

## Preface: Set up file system and input files
```{bash}
project_path=
mkdir -p ${project_path}/data/sequences/ref/
mkdir -p ${project_path}/data/sequences/query/
mkdir -p ${project_path}/data/variants/
mkdir -p ${project_path}/data/chromsize/
mkdir -p ${project_path}/data/annotations/
mkdir -p ${project_path}/analyses/tree
```

Download the following files and add them to the following directories. All except query files available on Cyverse: \
`${project_path}/data/sequences/ref/` \
B73 v5 reference sequence: Zm-B73-REFERENCE-NAM-5.0.fa \
Repeat masked B73 v5 reference sequence: Zm-B73-REFERENCE-NAM-5.0_rm.fa \
`${project_path}/data/sequences/query/` \
From Phytozome: \
Acomosus: Acomosus_321_v3.softmasked.fa.gz \
Osativa: Osativa_323_v7.0.softmasked.fa.gz \
Taestivum: Taestivum_296_v2.softmasked.fa.gz \
Bdistachyon: Bdistachyon_556_v3.0.softmasked.fa.gz \
Sbicolor: Sbicolor_313_v3.0.softmasked.fa.gz \
Othomaeum: Othomaeum_386_v1.0.fa.gz \
Stuberosum: Stuberosum_448_v4.03.fa.gz \
Vvinifera: Vvinifera_457_Genoscope.12X.fa.gz \
From Ensembl: \
Athaliana: Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa.gz \
We used four additional genomes which are currently unreleased. Chrysopogon serrulatus and Hyparrhenia diplandra (Song et al. 2020 preprint), and genomes for Coelorachis and Vossia.\
`${project_path}/data/variants/` \
NAM founders SVs in a bed file: NAM_founders_SVs.bed \
NAM founders SNPs in a bed file: B73v5.NAM-illumina-snps-only_filtered.frq.bed \
`${project_path}/data/chromsize/` \
Size of all chromosomes in B73: Zm-B73-REFERENCE-NAM-5.0.size.bed \
`${project_path}/data/annotations/` \
B73 v5 genes: zea_maysb73_core_3_87_1.gff \
Open chromatin from merging leaf and inflorescence data from Ricci et al. 2019: open_chromatin.bed \
GFF3 of repeat masked bases in B73 v5: Zm-B73-REFERENCE-NAM-5.0-repeatmask.gff3.gz \
Genetic map from Ogut et al. 2015 lifted over to AGPv5: ogut_fifthcM_map_agpv5.bed\

Run all of the following scripts from within your `${project_path}`.\

## 1) Align genomes

### Step 1: Run lastdb.sh on the reference genome

```{bash, eval=FALSE}
bash ./scripts-alignment/lastdb.sh
```

```{bash}

ref=./data/sequences/ref/Zm-B73-REFERENCE-NAM-5.0.fa
ref_name=$( basename $ref )
db_prefix=./analyses/last/lastdb/${ref_name}-MAM4

mkdir -p ./analyses/last/lastdb/
lastdb -P0 -uMAM4 -R01 $db_prefix $ref
```

### Step 2: Make alignments

Run the following for each query genome in ${project_path}/data/sequences/query/
```{bash, eval = FALSE}
bash make_alignments.sh query
```
```{bash, eval=FALSE}

ref=./data/sequences/ref/Zm-B73-REFERENCE-NAM-5.0.fa
ref_name=$( basename $ref )
db_prefix=./analyses/last/lastdb/${ref_name}-MAM4
query=$1
query_name=$(basename $query .gz)
query_name=$(basename $query_name .fa)

mkdir -p ./analyses/last/mat/
mat=./analyses/last/mat/"$ref_name"_"$query_name".mat

last-train -P0 --revsym --matsym --gapsym -E0.05 -C2 $db_prefix $query > $mat

mkdir -p ./analyses/last/maf/

maf=./analyses/last/maf/Zm-B73-REFERENCE-NAM-5.0_"$query_name".maf

lastal -m50 -E0.05 -C2 -p $mat $db_prefix $query > $maf
# m50 makes allowed multiplicity of initial hits 50, each match lengthened
# until it occurs at most this many times
# E is threshold for errors
# C2 makes lastal discard any gapless alignment in 2 or more gapped alignments
# p specifies match/mismatch score matrix
# Then reference db
# And finally the query fasta


mkdir -p ./analyses/last/axt
axt=./analyses/last/axt/"$ref_name"_"$query_name".axt

maf-convert axt $maf > $axt

mkdir -p ./analyses/last/chain
chain=./analyses/last/chain/"$ref_name"_"$query_name".chain
merged_chain=./analyses/last/chain/"$ref_name"_"$query_name".all.chain

axtChain $axt $ref $query $chain -linearGap=loose -faQ -faT

# axtChain options
# linearGap=loose for distant species
# -faQ The specified qNibDir is a fasta file with multiple sequences for query
# -faT The specified tNibDir is a fasta file with multiple sequences for target

mkdir -p ./analyses/last/chain_merged
merged_chain=./analyses/last/chain/"$ref_name"_"$query_name".all.chain

chainMergeSort $chain > $merged_chain

mkdir -p ./analyses/chromsize/
ref_size=./analyses/chromsize/"$ref_name".size
query_size=./analyses/chromsize/"$query_name".size

if [ ! -f $ref_size ]; then
  faSize $ref -detailed > $ref_size
fi

if [ ! -f $query_size ]; then
  faSize $query -detailed > $query_size
fi

# note: if your chromosomes are id'd with numbers, faSize will add everything
# before .fa in the file name to the beginning of each chromosome id in the size
# file. This causes issues in the next step. Either make all id's non-numeric
# (e.g. 'chr1' instead of '1') or go in to size file and remove prefix in front
# of chromosome ids.

mkdir -p ./analyses/last/chain_prenet/
chain_prenet=./analyses/last/chain_prenet/"$ref_name"_"$query_name".all.pre.chain
chainPreNet $merged_chain $ref_size $query_size $chain_prenet

mkdir -p ./analyses/last/target_net/
mkdir -p ./analyses/last/query_net/
target_net=./analyses/last/target_net/"$ref_name"_"$query_name".net
query_net=./analyses/last/query_net/"$query_name"_"$ref_name".net
chainNet $chain_prenet $ref_size $query_size $target_net $query_net

# making 2bit files for netToAxt
mkdir -p ./analyses/last/2bit/
query_twobit=./analyses/last/2bit/"$query_name".2bit
ref_twobit=./analyses/last/2bit/"$ref_name".2bit

if [ ! -f $query_twobit ]; then
  faToTwoBit $query $query_twobit
fi

if [ ! -f $ref_twobit ]; then
  faToTwoBit $ref $ref_twobit
fi

mkdir -p ./analyses/last/net_axt/
net_axt=./analyses/last/net_axt/"$ref_name"_"$query_name".net.axt
netToAxt $target_net $chain_prenet $ref_twobit $query_twobit $net_axt

mkdir -p ./analyses/last/net_axt/net_maf
net_maf=./analyses/last/net_axt/net_maf/"$ref_name"_"$query_name".net.maf

axtToMaf $net_axt $ref_size $query_size $net_maf

# now to make mafs one to one

head -n 29 $maf > "$net_maf"_w_header
cat $net_maf >> "$net_maf"_w_header

one_maf=./analyses/last/net_axt/net_maf/"$ref_name"_"$query_name".1to1.maf

last-split -m1 "$net_maf"_w_header |
maf-swap |
awk -v q="$query_name" -v r="$ref_name" '/^s/ {$2 = (++s % 2 ? q "." : r ".") \
$2} 1' | last-split -m1 | \
maf-swap | last-postmask > $one_maf

```

### Step 3: Combine and filter maf files

```{bash, eval=FALSE}
bash ./scripts-alignment/process_mafs.sh
```

```{bash}

maf_array=($( ls -d ./analyses/last/net_axt/net_maf/*1to1.maf ))
combined_maf=./analyses/last/net_axt/net_maf/combined.maf

cat ${maf_array[@]:0:1} | sed -n '/##maf version=1 scoring=blastz/,$p' > \
${maf_array[@]:0:1}_tmp
cat ${maf_array[@]:1:1} | sed -n '/##maf version=1 scoring=blastz/,$p' > \
${maf_array[@]:1:1}_tmp

multiz ${maf_array[@]:0:1}_tmp ${maf_array[@]:1:1}_tmp 1 > $combined_maf

for maf in ${maf_array[@]:2};
do
  echo "processing " $maf
 cat $maf | sed -n '/##maf version=1 scoring=blastz/,$p' > \
 "$maf"_tmp
 multiz $combined_maf "$maf"_tmp 1 > "$combined_maf"_tmp
 mv "$combined_maf"_tmp $combined_maf
done

# and filter mafs so all blocks have Zea mays and are at least 20 bp long
mafFilter -minCol=20 -needComp="$ref_name" $combined_maf > "$combined_maf".filtered
```

### Step 4: Convert maf files to fasta files

```{bash}
bash ./scripts-alignment/maf_to_msa_fasta.sh
```

```{bash}

ref_rm=./data/sequences/ref/Zm-B73-REFERENCE-NAM-5.0_rm.fa

ref_rm_name=$( basename $ref_rm )

combined_maf=./analyses/last/net_axt/net_maf/combined.maf

# splitting combined maf by target sequence
mkdir -p ./analyses/last/split_maf/
outroot=./analyses/last/split_maf/
mafSplit -byTarget dummy.bed $outroot "$combined_maf".filtered -useFullSequenceName

# mafSplit throws an error if you don't put dummy.bed, even though byTarget
# means it's being ignored

path_to_phast=~/bin/phast
msa_view=${path_to_phast}/bin/msa_view
maf_dir=($( ls -d ./analyses/last/split_maf/*chr*.maf ))


mkdir -p ./analyses/last/msa_fasta/
mkdir -p ./data/sequences/ref/split/

faSplit byname $ref_rm ./data/sequences/ref/split/

path_to_match_masking=./scripts-alignment/matchMasking.pl

for maf_file in "${maf_dir[@]}"
do
  chr=$(basename $maf_file | sed -e 's/0\(.*\).maf/\1/')

  fasta=./analyses/last/msa_fasta/"$chr".fa
  ref_rm_chr=./data/sequences/ref/split/"$chr".fa
  rm_fasta=./analyses/last/msa_fasta/"$chr"_rm.fa
  $msa_view $maf_file -f -G 1 --refseq $ref_rm_chr > $fasta
  sed -E 's/> />/g' $fasta > "$fasta"_tmp && mv "$fasta"_tmp $fasta
  perl $path_to_match_masking \
  --ref $ref_rm_chr \
  --fasta $fasta \
  --out $rm_fasta
done
```

## 2) Identify conserved elements

### Step 5: Estimate neutral tree

```{bash}
project_path=
tree="((Acomosus_321_v3,((Osativa_323_v7,(Taestivum_296_v2,Bdistachyon_556_v3)),((Vossia,(Zm-B73-REFERENCE-NAM-5,(Coelorachis,(Chrysopogon_serrulatus_FLYE2,(Hyparrhenia_diplandra_flye2,Sbicolor_454_v3))))),Othomaeum_386_v1))),(Solanum_tuberosum,(Vitis_vinifera,Arabidopsis_thaliana)));"
Rscript ./scripts-gerp/estimate_neutral_tree.R ${project_path}/data/sequences/ref/split/ ${project_path}/analyses/last/msa_fasta/ ${project_path}/data/annotations/${gff_name} $tree ${project_path}/analyses/tree/neutral_tree.txt

```

```{r}
library("rphast")
args <- commandArgs(trailingOnly = TRUE)
ref_root <- args[1]
msa_root <- args[2]
feat_file <- args[3]
tree <- args[4]
tree_output <- args[5]

chr_vector <- paste("chr", c(1:10), sep = "")

chr <- chr_vector[1]

ref_file <- paste(ref_root, chr, ".fa", sep = "")

msa <- paste(msa_root, chr, ".fa", sep="")

align4d <- read.msa(msa, refseq = ref_file, format = "FASTA", do.4d = TRUE, features= read.feat(feat_file))

for (chr in chr_vector[2:length(chr_vector)]){
  ref_file <- paste(ref_root, chr, ".fa", sep = "")
  msa <- paste(msa_root, chr, ".fa", sep="")
  align4d2 <- read.msa(msa, refseq = ref_file, format = "FASTA", do.4d = TRUE, features= read.feat(feat_file))
  align4d <- concat.msa(list(align4d, align4d2))
}

neutralMod <- phyloFit(align4d, tree=tree, subst.mod="REV")
sink(tree_output)
neutralMod$tree
sink()
```

### Step 6: Run GERP

Run for all 10 maize chromosomes (so chr1, chr2, etc.)
```{bash}
bash ./scripts-gerp/run_gerp.sh chr
```

run_gerp.sh
```{bash}

chr=$1
mkdir -p {project_path}/analyses/gerp/
fasta=./analyses/last/msa_fasta/"$chr"_rm.fa


ref_rm=./data/sequences/ref/Zm-B73-REFERENCE-NAM-5.0_rm.fa
ref_rm_chr=./data/sequences/ref/split/"$chr".fa
gerpcol={path/to/GERPplusplus/gerpcol}
gerpelem={path/to/GERPplusplus/gerpelem}
tree=./analyses/tree/neutral_tree.txt
$gerpcol -f $rm_fasta -t $tree -v -e Zm-B73-REFERENCE-NAM-5 -j -a


mkdir -p ./analyses/gerp/
gerp=./analyses/gerp/"$chr"_rm.rates
mv "$rm_fasta".rates $gerp

$gerpelem -f $gerp

# # make bed file
#
awk -v chr="$chr" 'BEGIN {OFS="\t"}; {print chr, NR-1, NR, $1, $2}' $gerp > "$gerp".bed
# #
# # # bed file of positive scores
awk '$5 > 0' "$gerp".bed > "$gerp".pos.bed
awk -v chr="$chr" 'BEGIN {OFS="\t"} {print chr, $2-1, $3, $4, $5, $6, $7, $8}' \
"$gerp".elems > "$gerp".elems.bed
sort -k 1,1 -k2,2n "$gerp".elems.bed > "$gerp".elems.bed_tmp
mv "$gerp".elems.bed_tmp "$gerp".elems.bed
```

## 3) Statistical analyses

### Step 7: Enrichment analyses

```{bash}
bash ./scripts-enrichment/bed_enrichment.sh
```

bed_enrichment.sh
```{bash}

test_bed=./analyses/gerp/all.elems.bed
cat ./analyses/gerp/*.elems.bed >> $test_bed

SV=./data/variants/NAM_founders_SVs.bed
SNP=./data/variants/B73v5.NAM-illumina-snps-only_filtered.frq.bed
genome=./data/chromsize/Zm-B73-REFERENCE-NAM-5.0.size.bed

# working with v5 gff
# get genes
gff=./data/annotations/zea_maysb73_core_3_87_1.gff
sed -i -e 's/^/chr/' $gff
genes=./data/annotations/genes_B73v5.gff3
genes_merged=./data/annotations/genes_B73v5_merged.gff3

awk 'BEGIN {OFS="\t"} $3~/gene/ {print $0}' $gff |
awk 'BEGIN {OFS="\t"} $5>$4 {print $0}'> $genes
sort -k 1,1 -k4,4n $genes > "$genes"_tmp && mv "$genes"_tmp $genes
bedtools merge -i $genes > $genes_merged
cds=/home/aihudson/projects/nam/data/annotations/cds_B73v5.gff3
awk 'BEGIN {OFS="\t"} $3=="CDS" {print $0}' $gff |
awk 'BEGIN {OFS="\t"} $5>$4 {print $0}' > $cds

# coding - get bed file with just coding regions

gerp_cds=./analyses/gerp/gerp_cds.elems.bed
bedtools intersect -a $test_bed -b $cds > $gerp_cds
sort -k 1,1 -k2,2n $gerp_cds > "$gerp_cds"_tmp && mv "$gerp_cds"_tmp $gerp_cds
bedtools merge -i $gerp_cds > "$gerp_cds"_tmp && mv "$gerp_cds"_tmp $gerp_cds

# get non-genic regions
gerp_nongenic=./analyses/gerp/gerp_nongenic.elems.bed
bedtools intersect -v -a $test_bed -b $genes_merged > $gerp_nongenic
sort -k 1,1 -k2,2n $gerp_nongenic > "$gerp_nongenic"_tmp && mv "$gerp_nongenic"_tmp $gerp_nongenic
bedtools merge -i $gerp_nongenic > "$gerp_nongenic"_tmp && mv "$gerp_nongenic"_tmp $gerp_nongenic

# get non-coding regions
gerp_noncds=./analyses/gerp/gerp_noncds.elems.bed
bedtools intersect -v -a $test_bed -b $cds > $gerp_noncds
sort -k 1,1 -k2,2n $gerp_noncds > "$gerp_noncds"_tmp && mv "$gerp_noncds"_tmp $gerp_noncds
bedtools merge -i $gerp_noncds > "$gerp_noncds"_tmp && mv "$gerp_noncds"_tmp $gerp_noncds

sort -k 1,1 -k2,2n $SV > "$SV"_tmp && mv "$SV"_tmp $SV

SV_merged=./data/variants/NAM_founders_SVs_merged.bed
sort -k 1,1 -k2,2n $SNP > "$SNP"_tmp && mv "$SNP"_tmp $SNP
SNP_merged=./data/variants/B73v5.NAM-illumina-snps-only_filtered_merged.frq.bed
bedtools merge -i $SNP > $SNP_merged

sh_helper_subset=.scripts-enrichment/bed_enrichment_helper_subset.sh
sh_helper=./scripts-enrichment/bed_enrichment_helper.sh
r_helper=./scripts-enrichment/bed_enrichment_helper.R

# Are structural variants enriched in the bed file?
chi=./analyses/chi/gerp_sv_chi.txt
bash $sh_helper $test_bed $SV_merged $genome $chi
# #
# #
# #
Rscript $r_helper $chi

# Are deletions enriched in the bed file?
del=./data/variants/NAM_founders_del.bed
sort -k 1,1 -k2,2n $del > "$del"_tmp && mv "$del"_tmp $del
del_merged=./data/variants/NAM_founders_del_merged.bed
bedtools merge -i $del > $del_merged
chi=./analyses/chi/gerp_del_chi.txt
bash $sh_helper $test_bed $del_merged $genome $chi
Rscript $r_helper $chi

# Are insertions enriched in the bed file?
ins=./data/variants/NAM_founders_ins.bed
sort -k 1,1 -k2,2n $ins > "$ins"_tmp && mv "$ins"_tmp $ins
ins_merged=./data/variants/NAM_founders_ins_merged.bed
bedtools merge -i $ins > $ins_merged
chi=./analyses/chi/gerp_ne_ins_bp_chi.txt
bash $sh_helper $test_bed $ins_merged $genome $chi
Rscript $r_helper $chi


# Are translocation enriched in the bed file?
tra=./data/variants/NAM_founders_tra.bed
sort -k 1,1 -k2,2n $tra > "$tra"_tmp && mv "$tra"_tmp $tra
tra_merged=./data/variants/NAM_founders_tra_merged.bed
bedtools merge -i $tra > $tra_merged
chi=./analyses/chi/gerp_ne_tra_chi.txt
bash $sh_helper $test_bed $tra_merged $genome $chi
Rscript $r_helper $chi


# Are inversions enriched in the bed file?
inv=./data/variants/NAM_founders_inv.bed
sort -k 1,1 -k2,2n $inv > "$inv"_tmp && mv "$inv"_tmp $inv
inv_merged=./data/variants/NAM_founders_inv_merged.bed
bedtools merge -i $inv > $inv_merged
chi=./analyses/chi/gerp_inv_chi.txt
bash $sh_helper $test_bed $inv_merged $genome $chi
Rscript $r_helper $chi

# Are SNPs enriched in the bed file?

chi=./analyses/chi/gerp_snp_chi.txt
bash $sh_helper $test_bed $SNP $genome $chi
#
#
#
Rscript $r_helper $chi
#
# Are structural variants enriched in coding regions?

chi=./analyses/chi/gerp_cds_sv_chi.txt
gerp_cds=./analyses/gerp/gerp_cds.elems.bed
# bedtools intersect -a $test_bed -b $cds > $gerp_cds
# sort -k 1,1 -k2,2n $gerp_cds > "$gerp_cds"_tmp && mv "$gerp_cds"_tmp $gerp_cds
# bedtools merge -i $gerp_cds > "$gerp_cds"_tmp && mv "$gerp_cds"_tmp $gerp_cds
bash $sh_helper_subset $gerp_cds $SV_merged $genome $chi $test_bed
# #
# #
# #
Rscript $r_helper $chi


# # Are SNPs enriched in coding regions?

chi=./analyses/chi/gerp_cds_snp_chi.txt

bash $sh_helper_subset $gerp_cds $SNP $genome $chi $test_bed

Rscript $r_helper $chi

# Are SVs enriched in non-genic GERP regions?

chi=./analyses/chi/gerp_nongenic_sv_chi.txt

bash $sh_helper_subset $gerp_nongenic $SV_merged $genome $chi $test_bed

Rscript $r_helper $chi

# Are SNPs enriched in non-genic GERP regions?

chi=./analyses/chi/gerp_nongenic_snp_chi.txt
bash $sh_helper_subset $gerp_nongenic $SNP $genome $chi $test_bed

Rscript $r_helper $chi


# Are SVs enriched in non-coding regions?

chi=./analyses/chi/gerp_noncds_sv_chi.txt

bash $sh_helper_subset $gerp_noncds $SV $genome $chi $test_bed



Rscript $r_helper $chi

# Are SNPs enriched in non-coding regions?

chi=./analyses/chi/gerp_noncds_snp_chi.txt

bash $sh_helper_subset $gerp_noncds $SNP $genome $chi $test_bed
#
Rscript $r_helper $chi

mkdir -p ./analyses/windows
windows=./analyses/windows/windows
bedtools makewindows -b $genome -w 10000 > "$windows".bed

bedtools intersect -wao -a $windows -b $ins_merged > "$windows"_ins.bed
bedtools intersect -wao -a $windows -b $del_merged > "$windows"_del.bed
bedtools intersect -wao -a $windows -b $tra_merged > "$windows"_tra.bed
bedtools intersect -wao -a $windows -b $inv_merged > "$windows"_inv.bed
bedtools intersect -wao -a $windows -b $test_bed > "$windows"_gerp.bed
open=./data/annotations/open_chromatin.bed
bedtools intersect -wao -a $windows -b $open > "$windows"_open.bed
masked=./data/annotations/Zm-B73-REFERENCE-NAM-5.0-repeatmask.gff3.gz
masked_merged=./data/annotations/Zm-B73-REFERENCE-NAM-5.0-repeatmask_merged.bed
bedtools merge -i $masked > $masked_merged
bedtools intersect -wao -a $windows -b $masked > "$windows"_masked.bed


```

Helper scripts for Step 7:
bed_enrichment_helper.sh
```{bash}
test_bed_1=$1
test_bed_2=$2
genome=$3
output_file=$4

intersect_size=$(bedtools intersect -a $test_bed_1 -b $test_bed_2 | \
awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
test_bed_2_size=$(bedtools merge -i $test_bed_2 | \
awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
test_bed_1_size=$(bedtools merge -i $test_bed_1 | \
awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
genome_size=$(awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' $genome)
((test_bed_1_subtract_test_bed_2=$test_bed_1_size-$intersect_size))
((test_bed_2_subtract_test_bed_1=$test_bed_2_size-$intersect_size))
((genome_subtract_test_beds=$genome_size-$intersect_size-$test_bed_1_subtract_test_bed_2-$test_bed_2_subtract_test_bed_1))
echo -e "Intersect\t$intersect_size" > $output_file
echo -e "Test_bed_1_subtract_test_bed_2\t$test_bed_1_subtract_test_bed_2" >> $output_file
echo -e "Test_bed_2_subtract_test_bed_1\t$test_bed_2_subtract_test_bed_1" >> $output_file
echo -e "Genome_subtract_test_beds\t$genome_subtract_test_beds" >> $output_file

```

bed_enrichment_helper.R
```{r}
args <- commandArgs(trailingOnly = TRUE)
filename = args[1]

data <- read.table(filename)
a <- data$V2[1]
b <- data$V2[2]
c <- data$V2[3]
d <- data$V2[4]

results <- fisher.test(rbind(c(a,b),c(c,d)))

saveRDS(results, paste(filename, ".results.RDS", sep = ""))

sink(paste(filename, ".results", sep = ""))
print(results)
sink()

```

bed_enrichment_helper_subset.sh

```{bash}
test_bed_1=$1
test_bed_2=$2
genome=$3
output_file=$4
superset=$5

intersect_size=$(bedtools intersect -a $test_bed_1 -b $test_bed_2 | \
awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
test_bed_2_size=$(awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' $test_bed_2)
test_bed_1_size=$(awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' $test_bed_1)
genome_size=$(awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' $genome)
superset_test_bed_2_intersect=$(bedtools intersect -a $superset -b $test_bed_2 | \
awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
((test_bed_1_subtract_test_bed_2=$test_bed_1_size-$intersect_size))
((test_bed_2_subtract_test_bed_1=$test_bed_2_size-$superset_test_bed_2_intersect))
((genome_subtract_test_beds=$genome_size-$intersect_size-$test_bed_1_subtract_test_bed_2-$test_bed_2_subtract_test_bed_1))
echo -e "Intersect\t$intersect_size" > $output_file
echo -e "Test_bed_1_subtract_test_bed_2\t$test_bed_1_subtract_test_bed_2" >> $output_file
echo -e "Test_bed_2_subtract_superset\t$test_bed_2_subtract_test_bed_1" >> $output_file
echo -e "Genome_subtract_test_beds\t$genome_subtract_test_beds" >> $output_file

```


### Step 8: Regression analysis

```{bash}
Rscript ./scripts-enrichment/regression_analysis.R $project_path
```

regression_analysis.R
```{r}

args <- commandArgs(trailingOnly=TRUE)
project_path=args[1]
library(dplyr)
library(data.table)

# read in the genetic map
map <- read.table(file=paste(project_path, "data/annotations/ogut_fifthcM_map_agpv5.bed", sep=""))
colnames(map) <- c("chr", "start", "stop", "s_marker", "m_marker", "useless", "cm")
# read in windows with insertions
windows_ins <- fread(file=paste(project_path, "analyses/windows/windows_ins.bed", sep=""))
windows_ins <- windows_ins %>%
select(V1,V2,V3,V7)
colnames(windows_ins) <- c("chr", "start", "stop", "ins_overlap")
# get the number of insertions that overlap each window
windows_ins <- windows_ins %>%
  mutate(ins_number = ifelse(ins_overlap == 0, 0, 1)) %>%
  group_by(chr,start,stop) %>%
  summarise_each(list(sum=sum))

# reading in windows with GERP
windows_gerp <- fread(file=paste(project_path,"analyses/windows/windows_gerp.bed",sep=""))
windows_gerp <- windows_gerp %>%
select(V1,V2,V3,V12)
colnames(windows_gerp) <- c("chr", "start", "stop", "gerp_overlap")
# get overlap of each window with GERP elements
windows_gerp <- windows_gerp %>%
  group_by(chr,start,stop) %>%
  summarise_each(list(sum=sum))
colnames(windows_gerp)[4] <- "gerp_overlap"

windows <- merge(windows_ins,windows_gerp)

# read in windows with open chromatin
windows_open <- fread(file=paste(project_path,"analyses/windows/windows_open.bed",sep=""))
windows_open <- windows_open %>%
select(V1,V2,V3,V7)
colnames(windows_open) <- c("chr", "start", "stop", "open_overlap")
# get overlap of each window with open chromatin
windows_open <- windows_open %>%
  group_by(chr,start,stop) %>%
  summarise_each(list(sum=sum))
colnames(windows_open)[4] <- "open_overlap"

windows <- merge(windows,windows_open)

# read in windows with masked base pairs
windows_masked <- fread(file=paste(project_path, "analyses/windows/windows_masked.bed", sep=""))
windows_masked <- windows_masked %>%
select(V1,V2,V3,V7)
colnames(windows_masked) <- c("chr", "start", "stop", "masked_overlap")
# get overlap of each window with masked base pairs
windows_masked <- windows_masked %>%
  group_by(chr,start,stop) %>%
  summarise_each(list(sum=sum))
colnames(windows_masked)[4] <- "masked_overlap"

windows <- merge(windows,windows_masked)

windows <- windows %>%
group_by(chr)

map <- map %>%
group_by(chr)

# shift all genetic distances forward by size of most negative distance to make all distances positive
make_positive <- function(x, na.rm = FALSE) (x - min(x))
map <- map %>%
mutate_at(c("cm"), make_positive) %>%
filter(cm == cummax(cm)) %>%
filter(start == cummax(start))

windows_list <- list()

# impute the genetic distance of each window using spline command
for(i in unique(windows$chr)){
  print(i)
  windows_subset <- windows %>%
  filter(chr==i) %>%
  arrange(chr,start)
  print(dim(windows_subset))
  map_subset <- map %>%
  filter(chr==i) %>%
  arrange(chr,start)
  print(dim(map_subset))
  windows_subset$start_cm <-
  spline(x=map_subset$start,y=map_subset$cm,
    xout=windows_subset$start, method = "hyman")$y
  windows_subset$stop_cm <- spline(x=map_subset$start,y=map_subset$cm,
    xout=windows_subset$stop, method = "hyman")$y
  windows_subset <- windows_subset %>%
  mutate(cm=stop_cm-start_cm) %>%
  filter(cm<.15) %>%
  filter(cm>0)
  windows_list[[i]] <- windows_subset
}

windows <- do.call(rbind, windows_list)

# add deletions
windows_del <- fread(file=paste(project_path, "analyses/windows/windows_del.bed", sep=""))
windows_del <- windows_del %>%
select(V1,V2,V3,V7)
colnames(windows_del) <- c("chr", "start", "stop","del_overlap")
windows_del <- windows_del %>%
  mutate(del_number = ifelse(del_overlap == 0, 0, 1)) %>%
  group_by(chr,start,stop) %>%
  summarise_each(list(sum=sum))
windows <- merge(x = windows, y = windows_del, all.x = TRUE)

# add inversions
windows_inv <- fread(file=paste(project_path,"analyses/windows/windows_inv.bed", sep = ""))
windows_inv <- windows_inv %>%
select(V1,V2,V3,V7)
colnames(windows_inv) <- c("chr", "start", "stop","inv_overlap")
windows_inv <- windows_inv %>%
  mutate(inv_number = ifelse(inv_overlap == 0, 0, 1)) %>%
  group_by(chr,start,stop) %>%
  summarise_each(list(sum=sum))
windows <- merge(x = windows, y = windows_inv, all.x = TRUE)

saveRDS(windows, paste(project_path, "analyses/windows/windows_w_sv_number.RDS", sep=""))

# convert cm (which is in cm/kb) to cm/Mb
data <- data %>%
  mutate(cm_per_mb = cm*100)

# second part of figure: insertions, inversions, deletions, regressed on recombination rate
# make a deta.frame to hold the coefficients from the regression models
betas <- data.frame(SV=NA,Model=NA,Beta=NA, SE=NA)
newdata2 <- data.frame(
  cm_per_mb = rep(seq(from = min(data$cm_per_mb), to =
                           max(data$cm_per_mb), length.out = 100), 3),
  gerp_overlap = (rep(0, each = 100)),
  open_overlap = rep(0, each = 100),
  masked_overlap = rep(0, each = 100))
newdata2 <- newdata2 %>%
  mutate(length=10000)
# regression for insertions
m_ins <- glm(ins_number_sum ~ gerp_overlap + cm_per_mb + open_overlap + masked_overlap +
               open_overlap:gerp_overlap +
               open_overlap:cm_per_mb +
               gerp_overlap:cm_per_mb +
               masked_overlap:gerp_overlap + masked_overlap:cm_per_mb +
               masked_overlap, data = data, family="quasipoisson")
betas[1,] <- c("Insertions", "Full Model", as.numeric(m_ins$coefficients[2]), sqrt(vcov(m_ins))[2,2])

# predict number of insertions in window based on varying recombination rate
newdata2 <- cbind(newdata2, predict(m_ins, newdata2, type = "response", se.fit=TRUE))
newdata2 <- within(newdata2, {
  Insertion <- (fit)
  Insertion.LL <- (fit - 1.96 * se.fit)
  Insertion.UL <- (fit + 1.96 * se.fit)
})
newdata2 <- newdata2 %>%
  dplyr::select(-c(fit,se.fit,residual.scale))
# regression for inversions
m_inv <- glm(inv_number_sum ~ gerp_overlap + cm_per_mb + open_overlap + masked_overlap +
               open_overlap:gerp_overlap +
               open_overlap:cm_per_mb +
               gerp_overlap:cm_per_mb +
               masked_overlap:gerp_overlap + masked_overlap:cm_per_mb +
               masked_overlap, data = data, family="quasipoisson")
betas[3,] <- c("Inversions", "Full Model", as.numeric(m_inv$coefficients[2]), sqrt(vcov(m_inv))[2,2])
# prediction for inversions
newdata2 <- cbind(newdata2, predict(m_inv, newdata2, type = "response", se.fit=TRUE))
newdata2 <- within(newdata2, {
  Inversion <- (fit)
  Inversion.LL <- (fit - 1.96 * se.fit)
  Inversion.UL <- (fit + 1.96 * se.fit)
})
newdata2 <- newdata2 %>%
  dplyr::select(-c(fit,se.fit,residual.scale))
# regression for deletions
m_del <- glm(del_number_sum ~ gerp_overlap + cm_per_mb + open_overlap + masked_overlap +
               open_overlap:gerp_overlap +
               open_overlap:cm_per_mb +
               gerp_overlap:cm_per_mb +
               masked_overlap:gerp_overlap + masked_overlap:cm_per_mb +
               masked_overlap, data = data, family="quasipoisson")
betas[2,] <- c("Deletions", "Full Model", as.numeric(m_del$coefficients[2]), sqrt(vcov(m_del))[2,2])

# prediction for deletions
newdata2 <- cbind(newdata2, predict(m_del, newdata2, type = "response", se.fit=TRUE))
newdata2 <- within(newdata2, {
  Deletion <- (fit)
  Deletion.LL <- (fit - 1.96 * se.fit)
  Deletion.UL <- (fit + 1.96 * se.fit)
})

# making the figures

library(RColorBrewer)
library(tidyr)
library(ggplot2)

# multiplot code from: http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# plotting the coefficients of the regressions models
betas$Beta <- as.numeric(betas$Beta)
betas$SE <- as.numeric(betas$SE)
betas_filtered <- betas %>%
  filter(Model == "Full Model")
p1 <- ggplot(betas_filtered, aes(x=SV, y=Beta, colour=Model)) +
    geom_errorbar(aes(ymin=Beta-SE, ymax=Beta+SE), width=.1) +
    geom_point() +
  labs(tag="A") +
    geom_hline(yintercept = 0, linetype = "dashed")   +
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=22,face="bold"),
        legend.position = "none")

# Plotting the prediction for each type of SV by varying recombination rate
newdata2_ins <- newdata2 %>%
  dplyr::select(cm_per_mb,Insertion,Insertion.UL,Insertion.LL)
newdata2_ins$SV <- rep("Insertion", length(newdata2_ins$cm_per_mb))
newdata2_ins <- newdata2_ins %>%
  rename(Estimate = Insertion, UL = Insertion.UL, LL = Insertion.LL)
newdata2_del <- newdata2 %>%
  dplyr::select(cm_per_mb,Deletion,Deletion.UL,Deletion.LL)
newdata2_del$SV <- rep("Deletion", length(newdata2_del$cm_per_mb))
newdata2_del <- newdata2_del %>%
  rename(Estimate = Deletion, UL = Deletion.UL, LL = Deletion.LL)
newdata2_inv <- newdata2 %>%
  dplyr::select(cm_per_mb,Inversion,Inversion.UL,Inversion.LL)
newdata2_inv$SV <- rep("Inversion", length(newdata2_inv$cm_per_mb))
newdata2_inv <- newdata2_inv %>%
  rename(Estimate = Inversion, UL = Inversion.UL, LL = Inversion.LL)
newdata2_long <- rbind(newdata2_del,newdata2_inv,newdata2_ins)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

p2 <- ggplot(newdata2_long, aes(cm_per_mb, Estimate, fill = SV, color =SV)) +
  geom_line(size = .5) +
  labs(x = "cM per Mb", y = "Number overlapping", tag="B") +
  geom_ribbon(aes(ymin = LL, ymax = UL), alpha = .25)  +
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=22,face="bold")) +
  scale_fill_manual(values=cbPalette) +
  scale_color_manual(values=cbPalette) +
  theme(legend.title=element_text(size=14),
    legend.text=element_text(size=13))

multiplot(p1,p2)
```
