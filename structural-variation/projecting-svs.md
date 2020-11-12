# Projecting SVs to NAM lines

by Rafael Della Coletta and Candice Hirsch (September, 2019 - April, 2020)

> The goal of this analysis is to project structural variants (SVs) indentified in the NAM founders onto the RILs of each NAM population. To do this, we need both SNP and SV calls for the founders, and SNP data for all NAM lines.



<!-- TOC START min:1 max:3 link:true asterisk:false update:true -->
- [Projecting SVs to NAM lines](#projecting-svs-to-nam-lines)
  - [Project folder](#project-folder)
  - [Transfering data from CyVerse to local folder](#transfering-data-from-cyverse-to-local-folder)
  - [Requirements](#requirements)
  - [Resequencing data](#resequencing-data)
    - [Transform VCF into Hapmap format](#transform-vcf-into-hapmap-format)
    - [Collapse overlapping SVs](#collapse-overlapping-svs)
    - [Merge SV calls from SNIFFLES and BioNano](#merge-sv-calls-from-sniffles-and-bionano)
    - [Identifying SNPs that are within the boundaries of a SV](#identifying-snps-that-are-within-the-boundaries-of-a-sv)
  - [GBS data](#gbs-data)
    - [Creating hapmap files for each NAM population and removing SNPs within SVs](#creating-hapmap-files-for-each-nam-population-and-removing-snps-within-svs)
    - [Collapsing duplicated SNPs](#collapsing-duplicated-snps)
    - [Overlaying resequencing data into parental GBS data](#overlaying-resequencing-data-into-parental-gbs-data)
    - [Select best GBS markers](#select-best-gbs-markers)
    - [Summary](#summary)
  - [Merge SNPs with SVs](#merge-snps-with-svs)
  - [Projection](#projection)
  - [Merge all projected SVs of each population in one file](#merge-all-projected-svs-of-each-population-in-one-file)
  - [Upload final hapmap to Cyverse](#upload-final-hapmap-to-cyverse)
- [Projecting resequencing SNPs to NAM lines](#projecting-resequencing-snps-to-nam-lines)
  - [Remove parental SNPs within SVs for each family](#remove-parental-snps-within-svs-for-each-family)
  - [Prepare datasets for projection](#prepare-datasets-for-projection)
  - [Projection](#projection-1)
  - [Merge all projected SNPs of each family in one file](#merge-all-projected-snps-of-each-family-in-one-file)
  - [Upload final hapmap to Cyverse](#upload-final-hapmap-to-cyverse-1)
- [Creating SNP subsets for GWAS](#creating-snp-subsets-for-gwas)
  - [LD calculation](#ld-calculation)
  - [Subsets for GWAS](#subsets-for-gwas)
- [Full SNP and SV dataset](#full-snp-and-sv-dataset)
<!-- TOC END -->



## Project folder

All data, scripts, and output of analyses are located on the folder `/home/hirschc1/della028/projects/sv_nams/` from my account at the Minnesota Supercomputing Institute (MSI).

```bash
cd ~/projects/

mkdir -p sv_nams/{analysis,data,scripts}
```




## Transfering data from CyVerse to local folder

On November 27th, Dr. Arun Seetharam shared the data needed for SV projection via CyVerse (SNPs) and Slack (SVs). The SNP data was located in the CyVerse folder `/iplant/home/shared/NAM/PANDA/SVs-impute`. The following commands were used to transfer this data to my folder at MSI so I can do my analyses.


```bash
# go to data folder of the project
cd ~/projects/sv_nams/data

# log in to cyverse
iinit
# go to cyverse shared folder to download data
icd /iplant/home/shared/NAM/PANDA/SVs-impute
# check if files match what Arun described
ils
# download data
iget -K GBS-output.tar.gz
iget -K B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.vcf.gz

# decompress files
tar xvzf GBS-output.tar.gz
gunzip B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.vcf.gz
```

The SV dataset was sent via Slack, and they were two files from two different SV calling methods (SNIFFLES and BioNano), which were merged into a single file later. I downloaded the files on my personal Mac and transferred them to the `data` folder in my MSI account via FileZilla.

```bash
# decompress files
gunzip data/NAM-SVs_SNIFFLES-min100-max1Mb-depth25_v6.vcf.gz
gunzip data/NAM-SVs_BioNano-min1Mb.txt.gz
```

After downloading and decompressing the files, these are the data that I will be using:

* `NAM-SVs_SNIFFLES-min100-max1Mb-depth25_v6.vcf`: file with SNIFFLES SV calls (up to 1 Mb) for NAM founders.
* `NAM-SVs_BioNano-min1Mb.txt`: file with BioNano SV calls (curated SVs larger than 1 Mb) for NAM founders.
* `B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.vcf`: file with SNP calls for NAM founders.
* `GBS-output/populations.snps.vcf`: file with SNP calls (GBS) for all NAM lines.



## Requirements

| Software | Version | Additional libraries / modules                                                         |
| -------- | ------- | -------------------------------------------------------------------------------------- |
| R        | 3.3.3   | `data.table (v1.12.4)`, `ggplot2 (v3.2.0)`, `foreach (v1.4.7)`, `doParallel (v1.0.15)` |
| Python   | 3.6.6   | `argparse (v1.1)`, `pandas (v0.23.4)`, `natsort (v6.0.0)`                              |
| TASSEL   | 5.2.56  | -                                                                                      |
| vcftools | 0.1.17  | -                                                                                      |
| plink    | 1.9     | -                                                                                       |

> Note: most of the bash `for` loops below can all be parallelized for better perfomance using [GNU parallel](https://www.gnu.org/software/parallel/). I'm showing a sequential way of doing that because it's easier to understand.



## Resequencing data

### Transform VCF into Hapmap format

I transformed VCF files to hapmap files before any kind of analysis because they are much smaller, easier to parse and quicker to analyze. Additionally, hapmap format would be used for projections anyways. The software [TASSEL 5](https://www.maizegenetics.net/tassel) can do this relatively fast with SNPs. However, a small complication arises when doing that for structural variants, since hapmap files were originally designed to store genotypic information as nucleotides. The way we got around that was writing a custom script (`scripts/vcf2hapmap.py`) to consider each SV as binary data (i.e. either present or not present) and code them as "nucleotides". Thus, if a SV is `A`bsent in a genotype, it was coded as `AA`, but if the SV is `T`here, it was coded as `TT`.

Converting SNPs from original VCF to hapmap would a take lot of time using TASSEl because the file is ~23 Gb. Therefore, I performed a series of UNIX commands to split the VCF into chromosomes and scaffolds and transform each of them into hapmap format.


```bash
# go to project folder
cd ~/projects/sv_nams

# create folder to store intermediate files
mkdir data/tmp

# filter vcf files
for i in {1..10}; do
  grep "^#"  data/B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.vcf > data/tmp/NAM_founders_SNPs.chr$i.vcf
  grep -w "^chr$i" data/B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.vcf >> data/tmp/NAM_founders_SNPs.chr$i.vcf
done
# get scaffolds as well
grep "^#" data/B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.vcf > data/tmp/NAM_founders_SNPs.scaffs.vcf
grep "^scaf_" data/B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.vcf >> data/tmp/NAM_founders_SNPs.scaffs.vcf

# transform each filtered vcf in hapmap
for i in {1..10}; do
  run_pipeline.pl -Xmx10g -importGuess data/tmp/NAM_founders_SNPs.chr$i.vcf -export data/tmp/NAM_founders_SNPs.chr$i.hmp.txt -exportType HapmapDiploid
done

# do the same for scaffolds, but this needs to be sorted first (according to TASSEL)
run_pipeline.pl -SortGenotypeFilePlugin -inputFile data/tmp/NAM_founders_SNPs.scaffs.vcf -outputFile data/tmp/NAM_founders_SNPs.scaffs.sorted.vcf -fileType VCF
# transform to diploid format (e.g. "AA" instead of "A")
run_pipeline.pl -Xmx10g -importGuess data/tmp/NAM_founders_SNPs.scaffs.sorted.vcf -export data/tmp/NAM_founders_SNPs.scaffs.hmp.txt -exportType HapmapDiploid

# correct typo in a genotype: it's supposed to be M37W and not MS37W
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 scaffs; do
  sed -i 1s/MS37W/M37W/ data/tmp/NAM_founders_SNPs.$chr.hmp.txt
done
```

> TASSEL throws this error when sorting vcf files `ERROR net.maizegenetics.dna.map.PositionListBuilder - validateOrdering: Position	Chr:SCAF_100	Pos:79721	InsertionPos:0	Name:SSCAF_100_79721	Variants:A/C	MAF:NaN	Ref:A and Position	Chr:SCAF_99	Pos:109272	InsertionPos:0	Name:SSCAF_99_109272	Variants:T/A	MAF:NaN	Ref:T out of order`. However, I think it's just a warning showing which positions were in the wrong position. I'm able to load the sorted vcf file and transform it into hapmap format without problems. Also, no SNP is lost when sorting the file.

Then, I converted the two files with SV calls (SNIFFLES and BioNano) for all NAM founders into hapmap with two different python scripts, since each file has a different format.

```bash
# go to project folder
cd ~/projects/sv_nams

# for explanation on how to use the scripts...
python scripts/vcf2hapmap.py -h
python scripts/variants2hapmap.py -h

# convert SVs vcf to hmp
python scripts/vcf2hapmap.py data/NAM-SVs_SNIFFLES-min100-max1Mb-depth25_v6.vcf data/NAM_founders_SVs.sniffles.not-sorted.hmp.txt
python scripts/variants2hapmap.py data/NAM-SVs_BioNano-min1Mb.txt data/NAM_founders_SVs.bionano.not-sorted.hmp.txt

# sort hmp file
run_pipeline.pl -Xmx10g -SortGenotypeFilePlugin \
                -inputFile data/NAM_founders_SVs.sniffles.not-sorted.hmp.txt \
                -outputFile data/NAM_founders_SVs.sniffles.sorted.hmp.txt \
                -fileType Hapmap
run_pipeline.pl -Xmx10g -SortGenotypeFilePlugin \
                -inputFile data/NAM_founders_SVs.bionano.not-sorted.hmp.txt \
                -outputFile data/NAM_founders_SVs.bionano.sorted.hmp.txt \
                -fileType Hapmap
# convert to diploid format (need to exclude B73 from sniffles since it's not
# present in the bionano dataset)
run_pipeline.pl -Xmx10g -importGuess data/NAM_founders_SVs.sniffles.sorted.hmp.txt \
                -excludeTaxa B73 \
                -export data/NAM_founders_SVs.sniffles.hmp.txt,data/filter_json \
                -exportType HapmapDiploid
run_pipeline.pl -Xmx10g -importGuess data/NAM_founders_SVs.bionano.sorted.hmp.txt \
                -export data/NAM_founders_SVs.bionano.hmp.txt \
                -exportType HapmapDiploid
```

Additional information about the SV is displayed on its ID, since hapmap format doesn't have fields available for adding such information. For example, the ID `del.chr1.51711.71809` on the first column of the hapmap file means that the SV is a deletion on chr1 that starts at 51,711 and ends at 71,809. The second column will also have the chromosome location for that deletion, but the third column will contain the **midpoint position** for that SV (i.e. 61,760). These two columns will always be the coordinates according to the reference genome. Although there will be somewhat redundant information on IDs of most SVs (like DELs, DUPs, INSs, and INVs), the ID will contain very important info about translocations, as it will show the respective location of the TRA in the **non-reference chromosome**.

Importantly, any SV called as heterozygous in the VCF file (i.e. `0/1`) was considered as **not** having a SV, therefore they were coded as `AA`.



### Collapse overlapping SVs

One problem with the BioNano SV calls is that the same SV can have slightly different, but overlapping, boundaries. The boundaries are fuzzy because of the methodology itself that cannot accurately tell what the start and end position of the SV are. Thus, we decided to collapse the information of SVs that have overlapping boundaries as they may actually represent the same SV. The approach for collapsing was: for each SV type that had overlapping boundaries, I selected the smallest start position and the greatest end position to be the coordinates for the collapsed SV. Then, I just merged their genotypic calls. I wrote `scripts/collapse_bionano_SVs.R` to do this.

```bash
module load R
# go to project folder
cd ~/projects/sv_nams

# collapse SVs
Rscript scripts/collapse_bionano_SVs.R data/NAM_founders_SVs.bionano.hmp.txt
```

SV calls by SNIFFLES also had similar problem. The difference is that same genotype can have different calls for the overlapping SV (which was not the case for bionano; i.e. all bionano overlapping SVs didn't have conflicts within a parent). Thus, in order to identify and remove these duplicated SVs in the SNIFFLES calls, I wrote `scripts/collapse_sniffles_SVs.R`. The conditions for collapsing or not overlapping SVs are as follows:

* If calls for overlapping SVs **are the same** across all parents:
  - Keep only the SV with less missing data;
  - If they have the same amount of missing data; keep the largest SV.
* If calls for overlapping SVs **are not the same** across all parents:
  - Keep both SVs (i.e. don't collapse)

```bash
module load R
# go to project folder
cd ~/projects/sv_nams/

Rscript scripts/collapse_sniffles_SVs.R data/NAM_founders_SVs.sniffles.hmp.txt
```

> Translocations were not considered in this filtering.




### Merge SV calls from SNIFFLES and BioNano

After collapsing SVs, I merged information from SNIFFLES and BioNano in a single file called `data/NAM_founders_SVs.hmp.txt`.

```bash
# go to project folder
cd ~/projects/sv_nams

# make sure the order of NAMs in hapmaps are the same
head -n 1 data/NAM_founders_SVs.sniffles.collapsed.hmp.txt data/NAM_founders_SVs.bionano.collapsed.hmp.txt
# they are

# merge files (use "sed 1d" to skip header of the bionano hapmap)
cat data/NAM_founders_SVs.sniffles.collapsed.hmp.txt > data/NAM_founders_SVs.not-sorted.hmp.txt
sed 1d data/NAM_founders_SVs.bionano.collapsed.hmp.txt >> data/NAM_founders_SVs.not-sorted.hmp.txt
# sort hmp file
run_pipeline.pl -Xmx10g -SortGenotypeFilePlugin \
                -inputFile data/NAM_founders_SVs.not-sorted.hmp.txt \
                -outputFile data/NAM_founders_SVs.sorted.hmp.txt \
                -fileType Hapmap
# convert to diploid format
run_pipeline.pl -Xmx10g -importGuess data/NAM_founders_SVs.sorted.hmp.txt \
                -export data/NAM_founders_SVs.hmp.txt \
                -exportType HapmapDiploid
```



### Identifying SNPs that are within the boundaries of a SV

SNPs that are found inside deletions are problematic, because they will have segregation issues when you compare multiple lines that have or not that SV. Thus, I wrote `scripts/generate_SV_bed.py` to create a BED file with start and end positions of deletions smaller than 100kb. I set up a 100 kb threshold because there were some extremely large deletions (>100 Mb) that would make me remove nearly all SNPs in this step and more than ~95% of the deletions were within that range. Importantly, each NAM population will be filtered separately.

> Translocations can cause SNP segregation issues as well. However, dealing with translocations is even more complicated, especially for SV projection and downstream GWAS, and we will ignore them here.

```bash
# go to project folder
cd ~/projects/sv_nams

# correct name of a NAM parent in the merged file
sed -i "s/IL14/IL14H/" data/NAM_founders_SVs.hmp.txt
# column numbers corresponding to NAM parents range from 13 to 37 (skip B73_Ab10)
for i in {13..37}; do
  # get NAM name
  NAM=$(head -n 1 data/NAM_founders_SVs.hmp.txt | cut -f $i)
  # create SV files for each NAM cross
  cut -f 1-12,$i data/NAM_founders_SVs.hmp.txt > data/NAM_founders_SVs_B73x$NAM.hmp.txt
  # create bed file with deletion boundaries for each population
  python scripts/generate_SV_bed.py data/NAM_founders_SVs_B73x$NAM.hmp.txt data/tmp/SNPs_to_remove_B73x$NAM.bed
  # # remove bed file's header for TASSEL compatibility
  # sed -i 1d data/tmp/SNPs_to_remove_B73x$NAM.bed
done

# after running the above, I noticed that 3 NAM parent names were slighlty different
# from the parental names in the GBS data (upper and lower case problme)
# so, I decided to change names of crosses now to avoid mismatches downstream
mv data/NAM_founders_SVs_B73xHP301.hmp.txt data/NAM_founders_SVs_B73xHp301.hmp.txt
mv data/NAM_founders_SVs_B73xIL14H.hmp.txt data/NAM_founders_SVs_B73xIl14H.hmp.txt
mv data/NAM_founders_SVs_B73xOh7b.hmp.txt data/NAM_founders_SVs_B73xOh7B.hmp.txt

mv data/tmp/SNPs_to_remove_B73xHP301.bed data/tmp/SNPs_to_remove_B73xHp301.bed
mv data/tmp/SNPs_to_remove_B73xIL14H.bed data/tmp/SNPs_to_remove_B73xIl14H.bed
mv data/tmp/SNPs_to_remove_B73xOh7b.bed data/tmp/SNPs_to_remove_B73xOh7B.bed
```


## GBS data

### Creating hapmap files for each NAM population and removing SNPs within SVs

The GBS file in VCF format contains information about millions of SNPs in all NAM lines (+5000), making it a very large file. Parsing it as it is would take way too much time. Thus, before doing any filtering, I split the VCF file by chromosome so I can run the filtering process for each chromosome in parallel.

```bash
# go to project folder
cd ~/projects/sv_nams

# create folder to save temporary files
mkdir data/GBS-output/tmp

# print header to each new chromosome vcf file
# awk will print only lines that start with # and quit after first mismatch
# (this is faster than using grep, which will go through the entire file and then quit)
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 scaffs; do
  awk '{if(/^#/)print;else exit}' data/GBS-output/populations.snps.vcf > data/GBS-output/tmp/NAM_rils_SNPs.$chr.vcf
done

# split vcf file by chromosome
module load parallel
for i in {1..10}; do
  sem -j +0 "grep -w '^chr$i' data/GBS-output/populations.snps.vcf >> data/GBS-output/tmp/NAM_rils_SNPs.chr$i.vcf"
done
sem --wait
grep "^scaf_" data/GBS-output/populations.snps.vcf >> data/GBS-output/tmp/NAM_rils_SNPs.scaffs.vcf
```

The bottleneck of performance now is that I have 5000+ lines in the vcf file, and transforming this file into hapmap format with TASSEL will take a lot of time. Thus I used [vcftools](https://vcftools.github.io/index.html) to create a vcf file for each NAM family (i.e. 25 vcf files with ~200 RILs each). But first, I had to create a file telling which RILs belong to which family based on the information on http://maizecoop.cropsci.uiuc.edu/nam-rils.php using `scripts/create_file_with_nam_rils_info.py`.

```bash
# go to project folder
cd ~/projects/sv_nams

# get list of all NAM RIL names (and parents) with GBS data
head -n 1 data/GBS-output/populations.sumstats.tsv | cut -f 2 > data/nam_ril_populations.txt
# rearrange information in a table
python scripts/create_file_with_nam_rils_info.py data/nam_ril_populations.txt

# filter vcf files for each chromosome...
for chr in {1..10}; do
  qsub -v CHR=$chr ~/projects/sv_nams/scripts/filter_nam_rils_and_snps_within_svs_vcf.sh
done
# ...and for scaffolds
{
  # skip header of "nam_ril_populations.txt" file
  read
  # set delimeter to tab
  IFS="\t"
  # read file line by line
  while read -r line; do
    # get name of the cross being parsed
    cross=$(echo $line | cut -f 1)
    # check if directory exists; if it doesnt, create one to store results
    [[ -d data/GBS-output/tmp/$cross ]] || mkdir -p data/GBS-output/tmp/$cross
    # transform line of the file into multiploe lines so that vcftools recognize 1 genotype to keep per line
    echo $line |  tr "\t" "\n" | tr "," "\n" > data/GBS-output/tmp/$cross/genotypes_to_keep.txt
    # use vcftools to filter a vcf file
    vcftools --vcf data/GBS-output/tmp/NAM_rils_SNPs.scaffs.vcf \
             --keep data/GBS-output/tmp/$cross/genotypes_to_keep.txt \
             --exclude-bed data/tmp/SNPs_to_remove_$cross.bed \
             --out data/GBS-output/tmp/$cross/NAM_rils_SNPs.$cross.scaffs.not-in-SVs \
             --recode \
             --recode-INFO-all
  done
} < "data/nam_ril_populations.txt"

# once the above is done, i have to sort each vcf file and export to hapmap format
cd ~/projects/sv_nams/data/GBS-output/tmp/

# commands for sorting and transforming to hapmap
for cross in $(ls -d B73x*); do
  qsub -v CROSS=$cross ~/projects/sv_nams/scripts/vcf2hmp_nam-rils_snps-not-in-svs.sh
done
```


### Collapsing duplicated SNPs

After taking a closer look at the hapmap files generated, I noticed that about many SNP positions were duplicated, i.e. there were multiple calls for the same SNP as seen in this example:

| chr | pos | B73 | Parent 2 | RIL 1 | RIL 2 | RIL 3 |
| --- | --- | --- | -------- | ----- | ----- | ----- |
| 1   | 100 | TT  | NN       | NN    | NN    | TA    |
| 1   | 100 | TT  | NN       | NN    | AA    | AA    |
| 1   | 100 | NN  | AA       | NN    | NN    | NN    |

To correct that, I collapsed these duplicates with `scripts/collapse_GBS_markers.R` by letting only the unique call for each RIL and converting to `NN` in case there was a among the calls for that RIL. The duplicated SNPs showed in the previous table would be collapsed into:

| chr | pos | B73 | Parent 2 | RIL 1 | RIL 2 | RIL 3 |
| --- | --- | --- | -------- | ----- | ----- | ----- |
| 1   | 100 | TT  | AA       | NN    | AA    | NN    |


```bash
cd ~/projects/sv_nams/data/GBS-output/tmp/

# collapse duplicated SNPs
for cross in $(ls -d B73x*); do
  for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 scaffs; do
  echo "Rscript ~/projects/sv_nams/scripts/collapse_GBS_markers.R $cross/NAM_rils_SNPs.$cross.$chr.not-in-SVs.hmp.txt $cross"
  done
done > ~/projects/sv_nams/scripts/commands_for_collapse-GBS-SNPs.txt

qsub ~/projects/sv_nams/scripts/collapse_GBS_markers.sh
```

After collapsing the duplicated SNPs, I merged the all hapmap files from each chromosome into one file.

```bash
cd ~/projects/sv_nams/data/GBS-output/tmp/

# merge hapmap files
for cross in $(ls -d B73x*); do
  cat $cross/NAM_rils_SNPs.$cross.chr1.not-in-SVs.collapsed.hmp.txt > $cross/NAM_rils_SNPs.$cross.not-in-SVs.not-imputed.hmp.txt
  for chr in chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 scaffs; do
    sed 1d $cross/NAM_rils_SNPs.$cross.$chr.not-in-SVs.collapsed.hmp.txt >> $cross/NAM_rils_SNPs.$cross.not-in-SVs.not-imputed.hmp.txt
  done
done

# check number of SNPs per pop
wc -l B73*/*not-in-SVs.not-imputed.hmp.txt
# ~2.7M (with slightly different number per population)
```



### Overlaying resequencing data into parental GBS data

The NAM parents were genotyped by both resequencing and GBS. Thus, it is possible that some SNP calls disagree between the two methods. To remove such SNPs, I overlayed the resequencing data into the GBS data and turned a SNP call into `NN` if there was a disagreement using `scripts/overlay_reseq-parental-SNPs_onto_GBS-data.R`. This script produces a hapmap for each cross containing **only the parental** data for that cross. Also, note that a number of SNPs in GBS data are not represented in the resequecing data, which can decrease even more the number of SNPs used for projection.

```bash
cd ~/projects/sv_nams/data/GBS-output/tmp/

for cross in $(ls -d B73x*); do
  echo "Rscript ~/projects/sv_nams/scripts/overlay_reseq-parental-SNPs_onto_GBS-data.R $cross/NAM_rils_SNPs.$cross.not-in-SVs.not-imputed.hmp.txt ~/projects/sv_nams/data/tmp/NAM_founders_SNPs.chr1.hmp.txt $cross $cross/NAM_gbs-parents_SNPs.$cross.not-in-SVs.reseq-overlay.hmp.txt"
done > ~/projects/sv_nams/scripts/commands_for_overlay-parental-SNPs.txt

qsub ~/projects/sv_nams/scripts/overlay_reseq-parental-SNPs_onto_GBS-data.sh

# check number of SNPs per pop
wc -l B73*/*not-in-SVs.reseq-overlay.hmp.txt
# ~1M (with slightly different number per population)
```

> Note: Apparently, the parent Tzi8 doesn't have gbs data. Therefore, I used the entire resequencing data for that parent.



### Select best GBS markers

GBS data contain a lot of missing data and also a lot of redundant information (many SNPs tightly linked). Besides increasing computation time with such big dataset, I also found in my preliminary analysis that using this raw GBS data had a strong negative impact on projections. Thus, we decided to filter this dataset by selecting only polymorphic SNPs, SNPs present in at least 30% of RILs, using a sliding window approach to remove incorrect calls (see [Huang et al, Genome Research, 2009](https://genome.cshlp.org/content/19/6/1068.abstract)) and removing SNPs with allele frequency < 0.4 or > 0.6. To do that, I ran `scripts/select_best_SNPs_per_pop.R`.

```bash
# go to data folder
cd ~/projects/sv_nams/data/GBS-output/tmp/

# # filter SNPs
# for cross in $(ls -d B73x*); do
#   Rscript ~/projects/sv_nams/scripts/select_best_SNPs_per_pop.R $cross $cross/NAM_rils_SNPs.$cross.not-in-SVs.not-imputed.hmp.txt $cross/NAM_gbs-parents_SNPs.$cross.not-in-SVs.reseq-overlay.hmp.txt ~/projects/sv_nams/analysis/qc/filter_best_SNPs --max_missing=0.3 --window_size=15 --window_step=1 --min_snps_per_window=5
# done

for cross in $(ls -d B73x*); do
  qsub -v CROSS=$cross ~/projects/sv_nams/scripts/select_best_SNPs_per_pop.sh
done

# check how many SNPs remained
wc -l B73x*/*.not-imputed.best-markers.hmp.txt
# range from 13k to 50k
```

**In summary**, from initial ~2.7 million SNPs for each population (after collapsing duplicated markers), ~1 million SNPs were present in the parental resequecing data, and the number of best SNPs selected varied from ~13k to ~50k depending on the population (mean of ~32k SNPs and median of ~31k).

> Note: The reason why there is `not-imputed` in the filename is because during preliminary tests I tried to impute GBS SNPs using FSFHap from TASSEL to decrease the number of missing data. But after selecting the best markers, I found the imputing SNPs didn't reduce much the missing data and it was actually causing some troubles later during SV projection. So we decided not to impute SNPs at this stage.



### Summary

**Before filtering**

I wrote `scripts/summary_raw_gbs.R` to plot some basic statistics of the raw GBS SNPs (after collapsing duplicates) such as total number of RILs, total number of SNPs, percentage of missing data and percentage of polymorphic SNPs for each population. It also plots the distribution of missing data per population, missing data per RIL, and missing data per SNP.

| Average number SNPs | Average missing data | Average polymorphic |
| ------------------- | -------------------- | ------------------- |
| 2,733,429           | 85.48%               | 29.32%              |

```bash
cd ~/projects/sv_nams/

# create new folder to store qc results
mkdir -p analysis/qc

# summarize data
Rscript scripts/summary_raw_gbs.R data/GBS-output/tmp analysis/qc/raw_gbs
```

In order to visualize how the markers are distributed along the chromosomes, I ploted karyotypes of 3 random RILs for each population with `scripts/plot_ril_karyotypes.R`. I used the chromosome coordinates from `analysis/qc/B73_RefGen_V4_chrm_info.txt` and the centromere positions from `analysis/qc/centromeres_Schneider-2016-pnas_v4.bed`, which are from B73 v4 reference genome.


```bash
cd ~/projects/sv_nams/data/GBS-output/tmp/

# create karyotypes
for cross in $(ls -d B73x*); do
  Rscript ~/projects/sv_nams/scripts/plot_ril_karyotypes.R ~/projects/sv_nams/analysis/qc/B73_RefGen_V4_chrm_info.txt ~/projects/sv_nams/analysis/qc/centromeres_Schneider-2016-pnas_v4.bed $cross ~/projects/sv_nams/analysis/qc/karyotypes/raw-gbs $cross/NAM_rils_SNPs.$cross.not-in-SVs.not-imputed.hmp.txt $cross/NAM_gbs-parents_SNPs.$cross.not-in-SVs.reseq-overlay.hmp.txt --rils=random --parents_in_data=TRUE --overlay_reseq=TRUE
done

# # run extra qc with TASSEL
# for cross in $(ls -d B73x*); do
#   run_pipeline.pl -Xmx6g -importGuess $cross/NAM_rils_SNPs.$cross.not-in-SVs.not-imputed.hmp.txt -GenotypeSummaryPlugin -endPlugin -export ~/projects/sv_nams/analysis/qc/raw_gbs/NAM_rils_SNPs_raw-gbs_$cross\_OverallSummary
#   (echo $cross && grep "Proportion Missing" ~/projects/sv_nams/analysis/qc/raw_gbs/NAM_rils_SNPs_raw-gbs_$cross\_OverallSummary1.txt) | tr "\n" "\t" | paste -s -d "\t" >> ~/projects/sv_nams/analysis/qc/raw_gbs/missing_data_raw-gbs.txt
# done
```


**After filtering**

When I filtered the raw gbs, `scripts/select_best_SNPs_per_pop.R` already produces some summary data. Thus, I just needed to generate the karyotypes for the same RILs used before filtering.

| Average number SNPs | Average missing data | Average polymorphic |
| ------------------- | -------------------- | ------------------- |
| 32,190              | 21.88%               | 100%                |

```bash
cd ~/projects/sv_nams/data/GBS-output/tmp/

# create karyotypes
for cross in $(ls -d B73x*); do
  # ugly way to get the names of rils used to plot karyotype before imputation
  rils=$(ls ~/projects/sv_nams/analysis/qc/karyotypes/raw-gbs/*$cross* | xargs -n 1 basename | cut -d "_" -f 2 | cut -d "." -f 1 | paste -s -d ",")
  # plot karyotypes for those rils
  Rscript ~/projects/sv_nams/scripts/plot_ril_karyotypes.R ~/projects/sv_nams/analysis/qc/B73_RefGen_V4_chrm_info.txt ~/projects/sv_nams/analysis/qc/centromeres_Schneider-2016-pnas_v4.bed $cross ~/projects/sv_nams/analysis/qc/karyotypes/best-markers $cross/NAM_rils_SNPs.$cross.not-in-SVs.not-imputed.best-markers.hmp.txt $cross/NAM_gbs-parents_SNPs.$cross.not-in-SVs.reseq-overlay.hmp.txt --rils=$rils --parents_in_data=TRUE --overlay_reseq=FALSE
done

# # run extra qc with TASSEL
# for cross in $(ls -d B73x*); do
#   run_pipeline.pl -Xmx6g -importGuess $cross/NAM_rils_SNPs.$cross.not-in-SVs.not-imputed.best-markers.hmp.txt -GenotypeSummaryPlugin -endPlugin -export ~/projects/sv_nams/analysis/qc/best-markers/NAM_rils_SNPs_best-markers_$cross\_OverallSummary
#   (echo $cross && grep "Proportion Missing" ~/projects/sv_nams/analysis/qc/best-markers/NAM_rils_SNPs_best-markers_$cross\_OverallSummary1.txt) | tr "\n" "\t" | paste -s -d "\t" >> ~/projects/sv_nams/analysis/qc/best-markers/missing_data_best-markers.txt
# done
```

Finally, I wrote `scripts/correct_SNP-names_rils.R` to make sure SNPs from RILs have the same name as the ones from parents. This script will generate the file `NAM_rils_SNPs.$cross.not-in-SVs.not-imputed.best-markers.correct-marker-names.hmp.txt`, which will be used later when merging SNPs with SVs.

```bash
cd ~/projects/sv_nams/data/GBS-output/tmp/

for cross in $(ls -d B73x*); do
  Rscript ~/projects/sv_nams/scripts/correct_SNP-names_rils.R $cross/NAM_gbs-parents_SNPs.$cross.not-in-SVs.reseq-overlay.hmp.txt $cross/NAM_rils_SNPs.$cross.not-in-SVs.not-imputed.best-markers.hmp.txt
done
```


## Merge SNPs with SVs

Now that I have selected which SNPs will be used to anchor projections, it's time to merge this SNP data with the SV calls into the same file for each population. I wrote `scripts/merge_SVs_and_SNPs.R` to do this task, and it generates a hapmap for the two parents of a cross (both SNP and SV calls come from resequecing), and another hapmap with the RILs of that cross (SNP calls come from GBS, and SV calls were set to `NN` since that's what will be projected later).

```bash
module load R/3.6.0

cd ~/projects/sv_nams/data/GBS-output/tmp/

# merge svs
for cross in $(ls -d B73x*); do
  Rscript ~/projects/sv_nams/scripts/merge_SVs_and_SNPs.R ~/projects/sv_nams/data/NAM_founders_SVs_$cross.hmp.txt $cross/NAM_gbs-parents_SNPs.$cross.not-in-SVs.reseq-overlay.hmp.txt $cross/NAM_rils_SNPs.$cross.not-in-SVs.not-imputed.best-markers.correct-marker-names.hmp.txt ~/projects/sv_nams/data/NAM_parents_SVs-SNPs.$cross.hmp.txt ~/projects/sv_nams/data/NAM_rils_SVs-SNPs.$cross.best-markers.not-projected.hmp.txt
done

# check that parents and rils have the same number of markers
for cross in $(ls -d B73x*); do
  wc -l ~/projects/sv_nams/data/NAM_parents_SVs-SNPs.$cross.hmp.txt
  wc -l ~/projects/sv_nams/data/NAM_rils_SVs-SNPs.$cross.best-markers.not-projected.hmp.txt
done
```

> Note: the B73 parent should not have any SV present, because all SVs were called against the B73 reference genome. However, since in practice there might be some miscalls, Arun used B73Ab10 as a negative control. However, this parent has a lot of missing data. This missing data can be problematic for projection. Since it's expected that B73Ab10 has very few SVs, I converted all `NN` into `AA` and left all non-missing calls as they were.




## Projection

All the projections of parental SVs into RILs will be done for each NAM population separately using the FILLIN plugin from TASSEL 5. In short, this program generate haplotypes for each parent of a cross (`-FILLINFindHaplotypesPlugin`), and then it imputes missing data in the RILs (i.e. SVs) using those parental haplotypes (`-FILLINImputationPlugin`).

But, before that, I need to make sure the files are sorted:

```bash
cd ~/projects/sv_nams/data/GBS-output/tmp/

# sort parents sv+snps
for cross in $(ls -d B73x*); do
  run_pipeline.pl -Xmx10g -SortGenotypeFilePlugin -inputFile ~/projects/sv_nams/data/NAM_parents_SVs-SNPs.$cross.hmp.txt -outputFile ~/projects/sv_nams/data/NAM_parents_SVs-SNPs.$cross.sorted.hmp.txt -fileType Hapmap
  run_pipeline.pl -Xmx10g -importGuess ~/projects/sv_nams/data/NAM_parents_SVs-SNPs.$cross.sorted.hmp.txt -export ~/projects/sv_nams/data/NAM_parents_SVs-SNPs.$cross.sorted.hmp.txt -exportType HapmapDiploid
done

# sort rils sv+snps
for cross in $(ls -d B73x*); do
  run_pipeline.pl -Xmx10g -SortGenotypeFilePlugin -inputFile ~/projects/sv_nams/data/NAM_rils_SVs-SNPs.$cross.best-markers.not-projected.hmp.txt -outputFile ~/projects/sv_nams/data/NAM_rils_SVs-SNPs.$cross.best-markers.not-projected.sorted.hmp.txt -fileType Hapmap
  run_pipeline.pl -Xmx10g -importGuess ~/projects/sv_nams/data/NAM_rils_SVs-SNPs.$cross.best-markers.not-projected.sorted.hmp.txt -export ~/projects/sv_nams/data/NAM_rils_SVs-SNPs.$cross.best-markers.not-projected.sorted.hmp.txt -exportType HapmapDiploid
done

# make sure number of rows of parental and RIL data matches in each population
for cross in $(ls -d B73x*); do
  wc -l ~/projects/sv_nams/data/NAM_parents_SVs-SNPs.$cross.sorted.hmp.txt
  wc -l ~/projects/sv_nams/data/NAM_rils_SVs-SNPs.$cross.best-markers.not-projected.sorted.hmp.txt
done
```

After some preliminary tests, I found that using a haplotype block size of 3,000 sites (`-hapSize 3000`) gives the best projection rate without compromising accuracy too much. Also, since I'm creating haplotypes for each parent individually, I need to set the minimum number of taxa to create a haplotype to 1 (`-minTaxa 1`). Importantly, the option `-hybNN` should be turned to `false`, otherwise the algorithm will combine haplotypes in recombination breakpoints if it thinks that region is heterozygous. If this option is not turned off, some projected regions will be messy and can end up with more than two alleles. It it's off, nothing is projected for that region.

```bash
# create new folder
mkdir ~/projects/sv_nams/analysis/projection

# project svs
for cross in $(ls -d B73x*); do
  qsub -v CROSS=$cross,HAPSIZE=3000 ~/projects/sv_nams/scripts/project_svs.sh
done
```

After projection, I wrote `scripts/count_projected_SVs.R` to compute how many SVs were projected per population and make summary plots about projections. I also wrote `scripts/plot_ril_karyotypes_SVs.R` to plot karyotypes of few RILs showing from which parent the SVs come from. The RILs selected were the same ones previously used to plot karyotypes of raw GBS data and best selected SNP markers.

```bash
cd ~/projects/sv_nams/data/GBS-output/tmp/

# calculate amount of projected SVs
qsub ~/projects/sv_nams/scripts/count_projected_svs.sh

# plot karyotypes of SVs present in each parent of a cross
for cross in $(ls -d B73x*); do
  Rscript ~/projects/sv_nams/scripts/plot_ril_karyotypes_SVs.R ~/projects/sv_nams/analysis/qc/B73_RefGen_V4_chrm_info.txt ~/projects/sv_nams/analysis/qc/centromeres_Schneider-2016-pnas_v4.bed $cross ~/projects/sv_nams/analysis/qc/karyotypes/projection ~/projects/sv_nams/analysis/projection/NAM_rils_SVs-SNPs.$cross.best-markers.projected.hmp.txt ~/projects/sv_nams/data/NAM_parents_SVs-SNPs.$cross.sorted.hmp.txt --rils=$rils --expected-SVs=TRUE
done

# now plot karyotypes with projected SVs for few RILs of each cross
for cross in $(ls -d B73x*); do
  # ugly way to get the names of rils used to plot karyotype before imputation
  rils=$(ls ~/projects/sv_nams/analysis/qc/karyotypes/best-markers/*$cross* | xargs -n 1 basename | cut -d "_" -f 2 | cut -d "." -f 1 | paste -s -d ",")
  # plot karyotypes for those rils
  Rscript ~/projects/sv_nams/scripts/plot_ril_karyotypes_SVs.R ~/projects/sv_nams/analysis/qc/B73_RefGen_V4_chrm_info.txt ~/projects/sv_nams/analysis/qc/centromeres_Schneider-2016-pnas_v4.bed $cross ~/projects/sv_nams/analysis/qc/karyotypes/projection ~/projects/sv_nams/analysis/projection/NAM_rils_SVs-SNPs.$cross.best-markers.projected.hmp.txt ~/projects/sv_nams/data/NAM_parents_SVs-SNPs.$cross.sorted.hmp.txt --rils=$rils --expected-SVs=FALSE
done

# additional QC about missing data
for cross in $(ls -d B73x*); do
  run_pipeline.pl -Xmx6g -importGuess ~/projects/sv_nams/analysis/projection/NAM_rils_SVs-SNPs.$cross.best-markers.projected.hmp.txt -GenotypeSummaryPlugin -endPlugin -export ~/projects/sv_nams/analysis/projection/NAM_rils_SVs-SNPs_$cross\_OverallSummary
  (echo $cross && grep "Proportion Missing" ~/projects/sv_nams/analysis/projection/NAM_rils_SVs-SNPs_$cross\_OverallSummary1.txt) | tr "\n" "\t" | paste -s -d "\t" >> ~/projects/sv_nams/analysis/projection/missing_data_best-markers_after_SV-projection.txt
done
# average missing data 0.09
```

The average percentage of projected SVs across all populations was **86%** (~135k SVs) with average accuracy of **93%**. The amount projected was a bit higher (87%) if considering only polymorphic SVs between the two parents of a cross. Only five crosses had projection rate below 75% (B73xCML322, B73xCML333, B73xOh7B, B73xP39, and B73xTzi8).



## Merge all projected SVs of each population in one file

The final file that will be sent to Jianming's group for GWAS will be a hapmap file **only with SVs** with all RILs of all NAM populations. I wrote `scripts/merge_SVs_after_projection.R` to filter and merge all SVs in one file. After combining multiple files, I have to correct the `alleles` column on the final hapmap to make sure that all possible alleles are represented in that column.

```bash
cd ~/projects/sv_nams/analysis/projection

# merge hapmaps
qsub ~/projects/sv_nams/scripts/merge_SVs_after_projection.sh

# make sure final file is sorted and let TASSEL correct the alleles' column
run_pipeline.pl -Xmx10g \
                -SortGenotypeFilePlugin \
                -inputFile ~/projects/sv_nams/analysis/projection/NAM_rils_projected-SVs-only.all-RILs.hmp.txt \
                -outputFile ~/projects/sv_nams/analysis/projection/NAM_rils_projected-SVs-only.all-RILs.final.hmp.txt \
                -fileType Hapmap
run_pipeline.pl -Xmx10g \
                -importGuess ~/projects/sv_nams/analysis/projection/NAM_rils_projected-SVs-only.all-RILs.final.hmp.txt \
                -export ~/projects/sv_nams/analysis/projection/NAM_rils_projected-SVs-only.all-RILs.final.hmp.txt \
                -exportType HapmapDiploid
```

After that, I performed a quick QC to check the total amount of missing data in the projected RILs.

```bash
cd ~/projects/sv_nams/analysis/projection

# get summary
run_pipeline.pl -Xmx40g -importGuess NAM_rils_projected-SVs-only.all-RILs.final.hmp.txt \
                -GenotypeSummaryPlugin -endPlugin \
                -export tassel_summary_NAM_rils_projected_svs

# plot missing data
Rscript ~/projects/sv_nams/scripts/qc_tassel_summary.R tassel_summary_NAM_rils_projected_svs3.txt \
                                                       missing_data_RILs_after_projection.png
```


## Upload final hapmap to Cyverse

Lastly, I uploaded the final hapmap file `~/projects/sv_nams/analysis/projection/NAM_rils_projected-SVs-only.all-RILs.final.v7.hmp.txt` to the shared folder on Cyverse.

```bash
# go to data folder of the project
cd ~/projects/sv_nams/analysis/projection

# change name to avoid conflict with previous versions
mv NAM_rils_projected-SVs-only.all-RILs.final.hmp.txt NAM_rils_projected-SVs-only.all-RILs.final.v7.hmp.txt

# log in to cyverse
iinit
# go to cyverse shared folder to download data
icd /iplant/home/shared/NAM/Misc
# check if files match what Arun described
ils
# upload data
iput -K NAM_rils_projected-SVs-only.all-RILs.final.v7.hmp.txt
# exit iRods
iexit full
```



# Projecting resequencing SNPs to NAM lines

The NAM RILs have a high amount of missing SNPs because they were called from GBS data. Since the whole genome resequencing data for the NAM parents include ~27 million SNPs, we also decided to project these SNPs into RILs in a similar manner as described above for SVs.


## Remove parental SNPs within SVs for each family

I previously removed SNPs within SVs for the RIL data. Now I need to also remove SNPs within the boundaries of deletions up to 100kb in the parental resequencing data (`data/tmp/NAM_founders_SNPs.vcf`) separately for each NAM family, since these are the SNPs that will be projected.

```bash
# go to project folder
cd ~/projects/sv_nams

# correct typos and case of parental names
# first get the vcf header
grep "^#" data/B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.vcf > data/tmp/NAM_founders_SNPs.vcf
# add substitution expressions to sed command
for sub in HP301/Hp301 IL14H/Il14H MS37W/M37W Ms71/MS71 OH43/Oh43 OH7B/Oh7B TX303/Tx303 TZi8/Tzi8; do
  sed -i "\$s/$sub/" data/tmp/NAM_founders_SNPs.vcf
done
# add the rest of data
grep -v "^#" data/B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.vcf >> data/tmp/NAM_founders_SNPs.vcf

# # remove SNPs within SVs and divide file into NAM families
# {
#   # skip header of "nam_ril_populations.txt" file
#   read
#   # set delimeter to tab
#   IFS="\t"
#   # read file line by line
#   while read -r line; do
#     # get name of the cross being parsed
#     cross=$(echo $line | cut -f 1)
#     echo "$cross"
#     # check if directory exists; if it doesnt, create one to store results
#     [[ -d data/tmp/$cross ]] || mkdir -p data/tmp/$cross
#     # add the two parents of a cross in a file so that vcftools recognize 1 genotype to keep per line
#     echo "B73" > data/tmp/$cross/genotypes_to_keep.txt
#     echo $cross | cut -d "x" -f 2-3 >> data/tmp/$cross/genotypes_to_keep.txt
#     # use vcftools to filter a vcf file
#     vcftools --vcf data/tmp/NAM_founders_SNPs.vcf \
#              --keep data/tmp/$cross/genotypes_to_keep.txt \
#              --exclude-bed data/tmp/SNPs_to_remove_$cross.bed \
#              --out data/tmp/$cross/NAM_parents-reseq_SNPs.$cross.not-in-SVs \
#              --recode
#   done
# } < "data/nam_ril_populations.txt"

qsub scripts/filter_nam_parents_and_snps_within_svs_reseq_vcf.sh
```

Once the above job is done, I have to sort each vcf file and export into diploid hapmap format:

```bash
cd ~/projects/sv_nams/data/tmp/

# commands for sorting
# run_pipeline.pl -Xmx10g -SortGenotypeFilePlugin -inputFile B73xB97/NAM_parents-reseq_SNPs.B73xB97.not-in-SVs.recode.vcf -outputFile B73xB97/NAM_parents-reseq_SNPs.B73xB97.not-in-SVs.sorted.vcf -fileType VCF
for cross in $(ls -d B73x*); do
  echo "run_pipeline.pl -Xmx10g -SortGenotypeFilePlugin -inputFile $cross/NAM_parents-reseq_SNPs.$cross.not-in-SVs.recode.vcf -outputFile $cross/NAM_parents-reseq_SNPs.$cross.not-in-SVs.sorted.vcf -fileType VCF"
done > ~/projects/sv_nams/scripts/commands_sort_vcf_reseq.txt

# commands for transforming to hapmap
# run_pipeline.pl -Xmx10g -importGuess B73xB97/NAM_parents-reseq_SNPs.B73xB97.not-in-SVs.sorted.vcf -export B73xB97/NAM_parents-reseq_SNPs.B73xB97.not-in-SVs.hmp.txt -exportType HapmapDiploid
for cross in $(ls -d B73x*); do
  echo "run_pipeline.pl -Xmx10g -importGuess $cross/NAM_parents-reseq_SNPs.$cross.not-in-SVs.sorted.vcf -export $cross/NAM_parents-reseq_SNPs.$cross.not-in-SVs.hmp.txt -exportType HapmapDiploid"
done > ~/projects/sv_nams/scripts/commands_vcf2hmp_reseq.txt

module load parallel
parallel --jobs 10 < ~/projects/sv_nams/scripts/commands_sort_vcf_reseq.txt
parallel --jobs 10 < ~/projects/sv_nams/scripts/commands_vcf2hmp_reseq.txt
```



## Prepare datasets for projection

Before projections, I will keep only the polymorphic parental SNPs for each NAM family to reduce computational time during projections (since they are the only informative SNPs anyways). The donor datasets will have the anchor markers and the polymorphic resequencing SNPs from the parents, while the RIL datasets will have the same anchor markers, but all polymorphic resequencing SNPs from parents will be set to missing data (`NN`).


```bash
# get only polymorphic snps
cd ~/projects/sv_nams/data/tmp/
for cross in $(ls -d B73x*); do
  qsub -v CROSS=$cross ~/projects/sv_nams/scripts/keep_poly_reseq-snps_only.sh
done

# merge svs
cd ~/projects/sv_nams/data/tmp/
for cross in $(ls -d B73x*); do
  qsub -v CROSS=$cross ~/projects/sv_nams/scripts/merge_SNPs-reseq_and_SNPs-SVs.sh
done
```

After running the above commands, I talked to Candy and we decided to use only SNPs (not SNPs and SVs) as anchors for projection. So I need to remove the SVs for each of the files created above.

```bash
# parents
for file in NAM_parents_SNPs-reseq_and_SVs-SNPs.*.poly.hmp.txt; do
  grep -v -P "^del|^dup|^ins|^inv|^tra" $file > ${file/SNPs-reseq_and_SVs-SNPs/SNPs-reseq_and_best-SNPs}
done

# rils
for file in B73x*/NAM_rils_SNPs-reseq_and_SVs-SNPs.*.poly.chr-*.not-projected.hmp.txt; do
  grep -v -P "^del|^dup|^ins|^inv|^tra" $file > ${file/SNPs-reseq_and_SVs-SNPs/SNPs-reseq_and_best-SNPs}
done
```

In order to speed up analysis, I will break the datasets into chromosomes so I can run more things in parallel:

```bash
# break parental data into chromosomes
for cross in $(ls -d B73x*); do
  for chr in {1..10}; do
    head -n 1 ~/projects/sv_nams/data/tmp/NAM_parents_SNPs-reseq_and_best-SNPs.$cross.poly.hmp.txt > ~/projects/sv_nams/data/tmp/NAM_parents_SNPs-reseq_and_best-SNPs.$cross.poly.chr-$chr.hmp.txt
    awk -v chr="$chr" '$3 == chr' ~/projects/sv_nams/data/tmp/NAM_parents_SNPs-reseq_and_best-SNPs.$cross.poly.hmp.txt >> ~/projects/sv_nams/data/tmp/NAM_parents_SNPs-reseq_and_best-SNPs.$cross.poly.chr-$chr.hmp.txt
  done
done

for cross in $(ls -d B73x*); do
  head -n 1 ~/projects/sv_nams/data/tmp/NAM_parents_SNPs-reseq_and_best-SNPs.$cross.poly.hmp.txt > ~/projects/sv_nams/data/tmp/NAM_parents_SNPs-reseq_and_best-SNPs.$cross.poly.chr-scaffs.hmp.txt
  grep -P "\tSCAF_" ~/projects/sv_nams/data/tmp/NAM_parents_SNPs-reseq_and_best-SNPs.$cross.poly.hmp.txt >> ~/projects/sv_nams/data/tmp/NAM_parents_SNPs-reseq_and_best-SNPs.$cross.poly.chr-scaffs.hmp.txt
done
```



## Projection

The projections will also done with FILLIN plugin from TASSEL for each NAM family separately, but with a bigger haplotype size since there are much more markers to project now. After some preliminary tests, using `-hapSize 70000` yielded the best projection rate and accuracy.

```bash
cd ~/projects/sv_nams/data/tmp/

# sort parents sv+snps
for cross in $(ls -d B73x*); do
  for chr in 1 2 3 4 5 6 7 8 9 10 scaffs; do
    echo "run_pipeline.pl -Xmx5g -SortGenotypeFilePlugin -inputFile ~/projects/sv_nams/data/tmp/NAM_parents_SNPs-reseq_and_best-SNPs.$cross.poly.chr-$chr.hmp.txt -outputFile ~/projects/sv_nams/data/tmp/NAM_parents_SNPs-reseq_and_best-SNPs.$cross.poly.chr-$chr.sorted.hmp.txt -fileType Hapmap"
  done
done > ~/projects/sv_nams/scripts/commands_sort-parents_reseq-SNPs.txt

# sort rils sv+snps
for cross in $(ls -d B73x*); do
  for chr in 1 2 3 4 5 6 7 8 9 10 scaffs; do
    echo "run_pipeline.pl -Xmx5g -SortGenotypeFilePlugin -inputFile ~/projects/sv_nams/data/tmp/$cross/NAM_rils_SNPs-reseq_and_best-SNPs.$cross.poly.chr-$chr.not-projected.hmp.txt -outputFile ~/projects/sv_nams/data/tmp/$cross/NAM_rils_SNPs-reseq_and_best-SNPs.$cross.poly.chr-$chr.not-projected.sorted.hmp.txt -fileType Hapmap"
  done
done > ~/projects/sv_nams/scripts/commands_sort-rils_reseq-SNPs.txt

# submit job
qsub ~/projects/sv_nams/scripts/sort_parents-rils_reseq_SNPs.sh

# make sure the number of SNPs match among datasets
for cross in $(ls -d B73x*); do
  for chr in 1 2 3 4 5 6 7 8 9 10 scaffs; do
    wc -l ~/projects/sv_nams/data/tmp/NAM_parents_SNPs-reseq_and_best-SNPs.$cross.poly.chr-$chr.sorted.hmp.txt
    wc -l ~/projects/sv_nams/data/tmp/$cross/NAM_rils_SNPs-reseq_and_best-SNPs.$cross.poly.chr-$chr.not-projected.sorted.hmp.txt
    echo ""
  done
done

# create folder to store results
mkdir ~/projects/sv_nams/analysis/reseq_snps_projection2/

# create haplotypes from parents
for cross in $(ls -d B73x*); do
  for chr in 1 2 3 4 5 6 7 8 9 10 scaffs; do
    echo "run_pipeline.pl -Xmx5g -FILLINFindHaplotypesPlugin -hmp ~/projects/sv_nams/data/tmp/NAM_parents_SNPs-reseq_and_best-SNPs.$cross.poly.chr-$chr.sorted.hmp.txt -o ~/projects/sv_nams/analysis/reseq_snps_projection2/donors_$cross-chr-$chr -hapSize 70000 -minTaxa 1"
  done
done > ~/projects/sv_nams/scripts/commands_donors_reseq-SNPs2.txt

# impute ril genotypes based on parental haplotypes
for cross in $(ls -d B73x*); do
  for chr in 1 2 3 4 5 6 7 8 9 10 scaffs; do
    echo "run_pipeline.pl -Xmx5g -FILLINImputationPlugin -hmp ~/projects/sv_nams/data/tmp/$cross/NAM_rils_SNPs-reseq_and_best-SNPs.$cross.poly.chr-$chr.not-projected.sorted.hmp.txt -d ~/projects/sv_nams/analysis/reseq_snps_projection2/donors_$cross-chr-$chr -o ~/projects/sv_nams/analysis/reseq_snps_projection2/NAM_rils_SNPs-reseq_and_best-SNPs.$cross.poly.chr-$chr.projected.hmp.txt -hapSize 70000 -accuracy -hybNN false"
  done
done > ~/projects/sv_nams/scripts/commands_project_reseq-SNPs2.txt

# submit job
qsub ~/projects/sv_nams/scripts/project_reseq_snps.sh
```

I will run the sliding window approach again (but with higher window size: 45-bp window, 1-bp step size, minimum of 15 markers per window) to correct possible errors around recombination breakpoints after projections. But first, I need to transform donor and projected datasets to hapmap diploid format.

```bash
# transform parents hapmap to diploid format
cd ~/projects/sv_nams/data/tmp/
for cross in $(ls -d B73x*); do
  for chr in {1..10}; do
    run_pipeline.pl -Xmx10g -importGuess NAM_parents_SNPs-reseq_and_best-SNPs.$cross.poly.chr-$chr.sorted.hmp.txt -export ~/projects/sv_nams/analysis/reseq_snps_projection2/NAM_parents_SNPs-reseq_and_best-SNPs.$cross.poly.chr-$chr.sorted.hmp.txt -exportType HapmapDiploid
  done
done

# transform rils hapmap to diploid
cd ~/projects/sv_nams/data/tmp/
for cross in $(ls -d B73x*); do
  for chr in {1..10}; do
    echo "run_pipeline.pl -Xmx10g -importGuess ~/projects/sv_nams/analysis/reseq_snps_projection2/NAM_rils_SNPs-reseq_and_best-SNPs.$cross.poly.chr-$chr.projected.hmp.txt -export ~/projects/sv_nams/analysis/reseq_snps_projection2/NAM_rils_SNPs-reseq_and_best-SNPs.$cross.poly.chr-$chr.projected.hmp.txt -exportType HapmapDiploid"
  done
done > ~/projects/sv_nams/scripts/commands_diploid-hmp_rils.txt

module load parallel
parallel --jobs 10 < ~/projects/sv_nams/scripts/commands_diploid-hmp_rils.txt

# submit jobs for sliding window
cd ~/projects/sv_nams/data/tmp/
for cross in $(ls -d B73x*); do
  for chr in {1..10}; do
    qsub -v CHR=$chr,CROSS=$cross ~/projects/sv_nams/scripts/sliding_window_reseq_snps.sh
  done
done

# merge chromosomes from projected hmp
for cross in $(ls -d B73x*); do
  cat ~/projects/sv_nams/analysis/reseq_snps_projection2/NAM_rils_SNPs-reseq_and_best-SNPs.$cross.poly.chr-1.projected.sliding-window.hmp.txt > ~/projects/sv_nams/analysis/reseq_snps_projection2/NAM_rils_SNPs-reseq_and_best-SNPs.$cross.poly.projected.hmp.txt
  for chr in 2 3 4 5 6 7 8 9 10 scaffs; do
    sed 1d ~/projects/sv_nams/analysis/reseq_snps_projection2/NAM_rils_SNPs-reseq_and_best-SNPs.$cross.poly.chr-$chr.projected.sliding-window.hmp.txt >> ~/projects/sv_nams/analysis/reseq_snps_projection2/NAM_rils_SNPs-reseq_and_best-SNPs.$cross.poly.projected.hmp.txt
  done
  echo "$cross done!"
done

# correct alleles' column of hapmap
for cross in $(ls -d B73x*); do
  run_pipeline.pl -Xmx10g -importGuess ~/projects/sv_nams/analysis/reseq_snps_projection2/NAM_rils_SNPs-reseq_and_best-SNPs.$cross.poly.projected.hmp.txt -export ~/projects/sv_nams/analysis/reseq_snps_projection2/NAM_rils_SNPs-reseq_and_best-SNPs.$cross.poly.projected.hmp.txt -exportType HapmapDiploid
done
```

After that, I can summarize projections of resequencing SNPs and plot karyotypes. The average percentage of projected SNPs across all populations was **92%%** with average accuracy of **97%**, and the karyotypes pretty much agree with those from the SV projections.

```bash
# calculate amount of projected SVs
qsub ~/projects/sv_nams/scripts/count_projected_reseq_snps.sh

# create folder to save karyotypes
mkdir -p ~/projects/sv_nams/analysis/qc/karyotypes/reseq_snps_projection2

# now plot karyotypes with projected SVs for few RILs of each cross
cd ~/projects/sv_nams/data/tmp/
for cross in $(ls -d B73x*); do
  # ugly way to get the names of rils used to plot karyotype before imputation
  rils=$(ls ~/projects/sv_nams/analysis/qc/karyotypes/best-markers/*$cross* | xargs -n 1 basename | cut -d "_" -f 2 | cut -d "." -f 1 | paste -s -d ",")
  # plot karyotypes for those rils
  Rscript ~/projects/sv_nams/scripts/plot_ril_karyotypes_reseq-SNPs.R \
          ~/projects/sv_nams/analysis/qc/B73_RefGen_V4_chrm_info.txt \
          ~/projects/sv_nams/analysis/qc/centromeres_Schneider-2016-pnas_v4.bed \
          $cross \
          ~/projects/sv_nams/analysis/qc/karyotypes/reseq_snps_projection2 \
          ~/projects/sv_nams/analysis/reseq_snps_projection2/NAM_rils_SNPs-reseq_and_best-SNPs.$cross.poly.projected.hmp.txt \
          ~/projects/sv_nams/data/tmp/$cross/NAM_parents-reseq_SNPs.$cross.poly.not-in-SVs.hmp.txt \
          --rils=$rils
done
```

Since the projections were performed only for polymorphic SNPs to reduce computational time, now I have to add back the monomorphic SNPs for each family.

```bash
cd ~/projects/sv_nams/data/tmp/

# create commands' file
for cross in $(ls -d B73x*); do
  for chr in {1..10}; do
    echo "Rscript ~/projects/sv_nams/scripts/add_mono-reseq-SNPs_after_projection.R $cross $chr ~/projects/sv_nams/data/tmp ~/projects/sv_nams/analysis/reseq_snps_projection2"
  done
done > ~/projects/sv_nams/scripts/commands_add_mono-reseq-SNPs2.txt

# submit job
qsub ~/projects/sv_nams/scripts/add_mono_reseq_snps.sh
```


## Merge all projected SNPs of each family in one file

Now I need add back snps that were not present in a certain cross to create a final file with information about all SNPs for all RILs (even if a SNP is missing in a RIL because it's not present in the parents).

```bash
# create empty hmp file with SNPS from all crosses for each chromosome
# (this will help merge snps from different crosses in the next script)
for chr in {1..10}; do
  qsub -v CHR=$chr ~/projects/sv_nams/scripts/create_empty_reseq-snps_file.sh
done

# add back snps not present in a cross
for chr in {1..10}; do
  qsub -v CHR=$chr ~/projects/sv_nams/scripts/merge_reseq-snps_per_chr.sh
done
```

Write final RIL file with all resequencing SNPs and SVs projected.

```bash
# merge all crosses with projected resequencing SNPs
for chr in {1..10}; do
  # exclude first columns for all crosses
  cd ~/projects/sv_nams/data/tmp/
  for cross in $(ls -d B73x*); do
    echo $chr $cross
    cut -f 1-11 --complement ~/projects/sv_nams/analysis/reseq_snps_projection2/NAM_rils_SNPs-reseq_and_best-SNPs.$cross.reseq-snps-all-crosses.chr-$chr.projected.hmp.txt > ~/projects/sv_nams/analysis/reseq_snps_projection2/NAM_rils_SNPs-reseq_and_best-SNPs.$cross.reseq-snps-all-crosses.chr-$chr.rils-only.projected.hmp.txt
  done
  # join all rils in one file (keep entire hmp file for cross B73xB97 though)
  cd ~/projects/sv_nams/analysis/reseq_snps_projection2
  paste NAM_rils_SNPs-reseq_and_best-SNPs.B73xB97.reseq-snps-all-crosses.chr-$chr.projected.hmp.txt \
        NAM_rils_SNPs-reseq_and_best-SNPs.B73xCML103.reseq-snps-all-crosses.chr-$chr.rils-only.projected.hmp.txt \
        NAM_rils_SNPs-reseq_and_best-SNPs.B73xCML228.reseq-snps-all-crosses.chr-$chr.rils-only.projected.hmp.txt \
        NAM_rils_SNPs-reseq_and_best-SNPs.B73xCML247.reseq-snps-all-crosses.chr-$chr.rils-only.projected.hmp.txt \
        NAM_rils_SNPs-reseq_and_best-SNPs.B73xCML277.reseq-snps-all-crosses.chr-$chr.rils-only.projected.hmp.txt \
        NAM_rils_SNPs-reseq_and_best-SNPs.B73xCML322.reseq-snps-all-crosses.chr-$chr.rils-only.projected.hmp.txt \
        NAM_rils_SNPs-reseq_and_best-SNPs.B73xCML333.reseq-snps-all-crosses.chr-$chr.rils-only.projected.hmp.txt \
        NAM_rils_SNPs-reseq_and_best-SNPs.B73xCML52.reseq-snps-all-crosses.chr-$chr.rils-only.projected.hmp.txt \
        NAM_rils_SNPs-reseq_and_best-SNPs.B73xCML69.reseq-snps-all-crosses.chr-$chr.rils-only.projected.hmp.txt \
        NAM_rils_SNPs-reseq_and_best-SNPs.B73xHp301.reseq-snps-all-crosses.chr-$chr.rils-only.projected.hmp.txt \
        NAM_rils_SNPs-reseq_and_best-SNPs.B73xIl14H.reseq-snps-all-crosses.chr-$chr.rils-only.projected.hmp.txt \
        NAM_rils_SNPs-reseq_and_best-SNPs.B73xKi11.reseq-snps-all-crosses.chr-$chr.rils-only.projected.hmp.txt \
        NAM_rils_SNPs-reseq_and_best-SNPs.B73xKi3.reseq-snps-all-crosses.chr-$chr.rils-only.projected.hmp.txt \
        NAM_rils_SNPs-reseq_and_best-SNPs.B73xKy21.reseq-snps-all-crosses.chr-$chr.rils-only.projected.hmp.txt \
        NAM_rils_SNPs-reseq_and_best-SNPs.B73xM162W.reseq-snps-all-crosses.chr-$chr.rils-only.projected.hmp.txt \
        NAM_rils_SNPs-reseq_and_best-SNPs.B73xM37W.reseq-snps-all-crosses.chr-$chr.rils-only.projected.hmp.txt \
        NAM_rils_SNPs-reseq_and_best-SNPs.B73xMo18W.reseq-snps-all-crosses.chr-$chr.rils-only.projected.hmp.txt \
        NAM_rils_SNPs-reseq_and_best-SNPs.B73xMS71.reseq-snps-all-crosses.chr-$chr.rils-only.projected.hmp.txt \
        NAM_rils_SNPs-reseq_and_best-SNPs.B73xNC350.reseq-snps-all-crosses.chr-$chr.rils-only.projected.hmp.txt \
        NAM_rils_SNPs-reseq_and_best-SNPs.B73xNC358.reseq-snps-all-crosses.chr-$chr.rils-only.projected.hmp.txt \
        NAM_rils_SNPs-reseq_and_best-SNPs.B73xOh43.reseq-snps-all-crosses.chr-$chr.rils-only.projected.hmp.txt \
        NAM_rils_SNPs-reseq_and_best-SNPs.B73xOh7B.reseq-snps-all-crosses.chr-$chr.rils-only.projected.hmp.txt \
        NAM_rils_SNPs-reseq_and_best-SNPs.B73xP39.reseq-snps-all-crosses.chr-$chr.rils-only.projected.hmp.txt \
        NAM_rils_SNPs-reseq_and_best-SNPs.B73xTx303.reseq-snps-all-crosses.chr-$chr.rils-only.projected.hmp.txt \
        NAM_rils_SNPs-reseq_and_best-SNPs.B73xTzi8.reseq-snps-all-crosses.chr-$chr.rils-only.projected.hmp.txt > NAM_rils_SNPs-reseq_and_best-SNPs.reseq-snps-all-crosses.chr-$chr.projected.hmp.txt
done

# merge projected svs with projected snps
# but first need to break parental data of SVs into chromosomes
for chr in {1..10}; do
  head -n 1 ~/projects/sv_nams/analysis/projection/NAM_rils_projected-SVs-only.all-RILs.final.v7.hmp.txt > ~/projects/sv_nams/analysis/reseq_snps_projection2/NAM_rils_projected-SVs-only.all-RILs.chr-$chr.hmp.txt
  awk -v chr="$chr" '$3 == chr' ~/projects/sv_nams/analysis/projection/NAM_rils_projected-SVs-only.all-RILs.final.v7.hmp.txt >> ~/projects/sv_nams/analysis/reseq_snps_projection2/NAM_rils_projected-SVs-only.all-RILs.chr-$chr.hmp.txt
done

# # make sure order of columns are the same
# head -n 1 ~/projects/sv_nams/analysis/reseq_snps_projection2/NAM_rils_projected-SVs-only.all-RILs.chr-1.hmp.txt > ~/projects/sv_nams/header_svs_dup-rem.txt
# head -n 1 ~/projects/sv_nams/analysis/reseq_snps_projection2/NAM_rils_SNPs-reseq_and_best-SNPs.reseq-snps-all-crosses.chr-1.projected.hmp.txt > ~/projects/sv_nams/header_snps2_dup-rem.txt
# diff ~/projects/sv_nams/header_svs_dup-rem.txt ~/projects/sv_nams/header_snps2_dup-rem.txt

# and then finally merge svs to snps
cd ~/projects/sv_nams/analysis/reseq_snps_projection2
for chr in {1..10}; do
  echo $chr
  sed 1d NAM_rils_projected-SVs-only.all-RILs.chr-$chr.hmp.txt >> NAM_rils_SNPs-reseq_and_best-SNPs.reseq-snps-all-crosses.chr-$chr.projected.hmp.txt
done

# sort snps + svs
cd ~/projects/sv_nams/analysis/reseq_snps_projection2
for chr in {1..10}; do
  qsub -v CHR=$chr ~/projects/sv_nams/scripts/sort_reseq-snps_and_svs.sh
done

# fix alleles column and transform to diploid hapmap
cd ~/projects/sv_nams/analysis/reseq_snps_projection2
for chr in {1..10}; do
  qsub -v CHR=$chr ~/projects/sv_nams/scripts/fix_allele_columns_reseq-snps.sh
done

# merge all chromosomes
cd ~/projects/sv_nams/analysis/reseq_snps_projection2
cp NAM_rils_SNPs-reseq_and_best-SNPs.reseq-snps-all-crosses.chr-1.projected.hmp.txt NAM_rils_SNPs-reseq_and_best-SNPs.projected.final.v7.hmp.txt

for chr in {2..10}; do
  echo $chr
  sed 1d NAM_rils_SNPs-reseq_and_best-SNPs.reseq-snps-all-crosses.chr-$chr.projected.hmp.txt >> NAM_rils_SNPs-reseq_and_best-SNPs.projected.final.v7.hmp.txt
done
```

## Upload final hapmap to Cyverse

Lastly, I uploaded the final hapmap file `~/projects/sv_nams/analysis/reseq_snps_projection2/NAM_rils_SNPs-reseq_and_SVs.projected.final.v7.hmp.txt` to the shared folder on Cyverse.

```bash
# go to data folder of the project
cd ~/projects/sv_nams/analysis/reseq_snps_projection2

# i realized that the I haven't renamed the file after adding the SVs
mv NAM_rils_SNPs-reseq_and_best-SNPs.projected.final.v7.hmp.txt NAM_rils_SNPs-reseq_and_SVs.projected.final.v7.hmp.txt
# compress for faster upload
gzip -c NAM_rils_SNPs-reseq_and_SVs.projected.final.v7.hmp.txt > NAM_rils_SNPs-reseq_and_SVs.projected.final.v7.hmp.txt.gz

# log in to cyverse
iinit
# go to cyverse shared folder to download data
icd /iplant/home/shared/NAM/Misc
# check if files match what Arun described
ils
# upload data
iput -K NAM_rils_SNPs-reseq_and_SVs.projected.final.v7.hmp.txt.gz
# exit iRods
iexit full
```



# Creating SNP subsets for GWAS

We want to understand whether or not SVs can capture phenotypic variation not explained by SNPs alone in GWAS, but since there are much more SNPs than SVs in our dataset, it would be an unfair comparison to run GWAS only with either type of marker. In addition, given the high density of SNPs, the level of LD between SNPs and SVs can lead to a biased result. Thus, we want to subsample SNPs from RILs according to their level of LD to SVs and also match the total number of SVs.


## LD calculation

First thing to do is to calculate LD between SNPs and SVs, because it will take some time given the amount of RILs and markers in the NAM dataset. I will use Plink v1.9 for that, and will only calculate LD between markers that are in a 100kb window. Then I will filter the output to have only LD between non-translocation SVs and SNPs.

```bash
# transform files from hapmap to plk format
cd ~/projects/sv_nams/analysis/reseq_snps_projection2
mkdir ld
for chr in {1..10}; do
  qsub -v CHR=$chr ~/projects/sv_nams/scripts/hmp2plk.sh
done

# calculate LD with ld-window 1M variants and ld-window-kb 100kb
for chr in {1..10}; do
  qsub -v CHR=$chr ~/projects/sv_nams/scripts/plink_ld_snps-svs_100kb-window.sh
done

# remove translocations, and keep only SNPs and SVs LD
cd /scratch.global/della028/hirsch_lab/ld_files
for chr in {1..10}; do
  qsub -v CHR=$chr ~/projects/sv_nams/scripts/plink_only_snp-sv-ld.sh
done
```

> For LD calculations, set up the `--ld-window` to the same number as your `--ld-window-kb`. The former option looks at a window depending on the number of variants, while the latter looks at physical distance (in kb). This is important, because say you have 1000 variants in a 10 kb window. If you set the options `--ld-window 100` and `--ld-window-kb 10`, it will stop calculating LD after it looks at the first 100 variants (so you don’t cover the whole 10kb window). Also `--make-founders` will calculate ld among all your lines in your dataset. Also set `--ld-window-r2 0` to make sure plink reports r2 value of zero or more (the default is 0.2).


## Subsets for GWAS

To reduce interference of missing data in downstream analyses, we will select only non-translocation SVs with more than 80% data across the families that had information in the founders, and SNPs that have more than 80% data across all families (because we want to make sure to use SNPs that have low amount of missing data, otherwise, the LD between an SV and a SNP can be biased).

```bash
cd ~/projects/sv_nams

# get name of all SVs after removing duplicates
cut -f 1 analysis/projection/NAM_rils_projected-SVs-only.all-RILs.final.v7.hmp.txt | sed 1d > analysis/projection/SV_names_after_removing_duplicates.txt
# get summary of the number of projected RILs per SV ("summary_projected_RILs_per_sv.txt")
Rscript scripts/count_projected_SVs.R ~/projects/sv_nams/data ~/projects/sv_nams/analysis/projection

# get names of all SNPs first
grep -v -P "^#" data/tmp/NAM_founders_SNPs.vcf | cut -f 1-2 | grep -v "scaf" | tr "\t" "_" | sed "s/^chr/S/" > data/tmp/all_SNP_names.txt
# get summary of the number of projected RILs per reseq SNP ("summary_projected_RILs_per_reseq-snp.txt")
Rscript scripts/count_projected_reseq-SNPs.R ~/projects/sv_nams/data/tmp ~/projects/sv_nams/analysis/reseq_snps_projection2

# keep only SNPs present in 80% or more families in the LD file
cd ~/projects/sv_nams/analysis/reseq_snps_projection2
# create directory to store file
mkdir -p ld/missing_data_filter
# run tassel to determine amount of missing data per marker
for chr in {1..10}; do
 qsub -v CHR=$chr ~/projects/sv_nams/scripts/tassel_summary_final_NAM_set.sh
done
# go to folder
cd ld/missing_data_filter/
# get number of column with info about missing data
head -n 1 tassel_summary_chr-10_NAM_rils_3.txt | tr '\t' '\n' | cat -n | grep "Proportion Missing"
# 33	Proportion Missing
# get marker names to keep
for chr in {1..10}; do
  sed 1d tassel_summary_chr-$chr\_NAM_rils_3.txt | awk '$33 <= 0.2' | cut -f 2 | grep -P "^S" > SNPs_low_missing-data_chr-$chr.txt
done
cat SNPs_low_missing-data_chr-1.txt > SNPs_low_missing-data.txt
for chr in {2..10}; do
  cat SNPs_low_missing-data_chr-$chr.txt >> SNPs_low_missing-data.txt
done

# filter by amount of missing data given there are information about the marker in the parents
cd ~/projects/sv_nams
Rscript scripts/filter_markers_by_missing_RILs.R analysis/projection/summary_projected_RILs_per_sv.txt \
                                                 analysis/projection/SV_names_after_removing_duplicates.txt \
                                                 analysis/reseq_snps_projection2/summary_projected_RILs_per_reseq-snp.txt \
                                                 analysis/reseq_snps_projection2/ld/missing_data_filter/SNPs_low_missing-data.txt \
                                                 0.8 \
                                                 analysis/reseq_snps_projection2/ld/missing_data_filter

# get names of all non-translocation SVs for each chromosome
cd ~/projects/sv_nams/analysis/reseq_snps_projection2
for chr in {1..10}; do
 grep -v -P "^S" ld/missing_data_filter/markers_to_keep_chr-$chr\_missing_filter.txt > ~/projects/sv_nams/data/subset-NAM-snps/SVs-to-keep.missing-filter.no-tra.chr$chr.txt
done
# merge all chromosomes
cd ~/projects/sv_nams/data/subset-NAM-snps/
cp SVs-to-keep.missing-filter.no-tra.chr1.txt SVs-to-keep.missing-filter.no-tra.txt
for chr in {2..10}; do
 cat SVs-to-keep.missing-filter.no-tra.chr$chr.txt >> SVs-to-keep.missing-filter.no-tra.txt
done

# get names of all SNPs with more than 80% data for each chromosome
cd ~/projects/sv_nams/analysis/reseq_snps_projection2
for chr in {1..10}; do
 grep -P "^S" ld/missing_data_filter/markers_to_keep_chr-$chr\_missing_filter.txt > ~/projects/sv_nams/data/subset-NAM-snps/SNPs-to-keep.missing-filter.chr$chr.txt
done
# merge all chromosomes
cd ~/projects/sv_nams/data/subset-NAM-snps/
cp SNPs-to-keep.missing-filter.chr1.txt SNPs-to-keep.missing-filter.txt
for chr in {2..10}; do
 cat SNPs-to-keep.missing-filter.chr$chr.txt >> SNPs-to-keep.missing-filter.txt
done
```

There were **169,793 non-translocation SVs** with more than 80% data among RILs that had information on the founders. Then, we randomly selected SNPs to match that number of filtered SVs (subset 1), selected a single SNP in highest LD for each SV (or closest one if more than one SNP with highest LD; subset 2), and randomly selected SNPs that were not in LD with any SV (R2 < 0.2; subset 3).

```bash
# there are 277,527 SVs in total to sample
wc -l ~/projects/sv_nams/data/subset-NAM-snps/SVs-to-keep.missing-filter.no-tra.chr*.txt
# 14424 /home/hirschc1/della028/projects/sv_nams/data/subset-NAM-snps/SVs-to-keep.missing-filter.no-tra.chr10.txt
# 27830 /home/hirschc1/della028/projects/sv_nams/data/subset-NAM-snps/SVs-to-keep.missing-filter.no-tra.chr1.txt
# 13358 /home/hirschc1/della028/projects/sv_nams/data/subset-NAM-snps/SVs-to-keep.missing-filter.no-tra.chr2.txt
# 18647 /home/hirschc1/della028/projects/sv_nams/data/subset-NAM-snps/SVs-to-keep.missing-filter.no-tra.chr3.txt
# 16229 /home/hirschc1/della028/projects/sv_nams/data/subset-NAM-snps/SVs-to-keep.missing-filter.no-tra.chr4.txt
# 18624 /home/hirschc1/della028/projects/sv_nams/data/subset-NAM-snps/SVs-to-keep.missing-filter.no-tra.chr5.txt
# 16421 /home/hirschc1/della028/projects/sv_nams/data/subset-NAM-snps/SVs-to-keep.missing-filter.no-tra.chr6.txt
# 13169 /home/hirschc1/della028/projects/sv_nams/data/subset-NAM-snps/SVs-to-keep.missing-filter.no-tra.chr7.txt
# 16273 /home/hirschc1/della028/projects/sv_nams/data/subset-NAM-snps/SVs-to-keep.missing-filter.no-tra.chr8.txt
# 14818 /home/hirschc1/della028/projects/sv_nams/data/subset-NAM-snps/SVs-to-keep.missing-filter.no-tra.chr9.txt

# make sure to use SNPs that have R2 calculated to an SV with more than 80% data
cd /scratch.global/della028/hirsch_lab/ld_files
for chr in {1..10}; do
  echo $chr
  sed 1d NAM_rils_SNPs-reseq_and_best-SNPs.reseq-snps-all-crosses.chr-$chr.projected.no-tra.snp-sv.ld | tr -s [:space:] | sed -e 's/^ //g' | tr " " "\t" | cut -f 3,6 | tr "\t" "\n" | grep -v -P "^del|^dup|^ins|^inv" | sort | uniq > SNPs_after_plink_ld.$chr.no-tra.snp-sv.txt
done

for chr in {1..10}; do
  echo $chr
  grep -Fxf SNPs_after_plink_ld.$chr.no-tra.snp-sv.txt ~/projects/sv_nams/data/subset-NAM-snps/SNPs-to-keep.missing-filter.chr$chr.txt > SNPs_after_plink_ld.$chr.no-tra.snp-sv.missing-filter.txt
done

# subsample SNPs by chromosome based on number of SVs above
shuf SNPs_after_plink_ld.1.no-tra.snp-sv.missing-filter.txt -n 27830 -o ~/projects/sv_nams/data/subset-NAM-snps/SNPs-to-keep_random-missing-filter_chr1.txt
shuf SNPs_after_plink_ld.2.no-tra.snp-sv.missing-filter.txt -n 13358 -o ~/projects/sv_nams/data/subset-NAM-snps/SNPs-to-keep_random-missing-filter_chr2.txt
shuf SNPs_after_plink_ld.3.no-tra.snp-sv.missing-filter.txt -n 18647 -o ~/projects/sv_nams/data/subset-NAM-snps/SNPs-to-keep_random-missing-filter_chr3.txt
shuf SNPs_after_plink_ld.4.no-tra.snp-sv.missing-filter.txt -n 16229 -o ~/projects/sv_nams/data/subset-NAM-snps/SNPs-to-keep_random-missing-filter_chr4.txt
shuf SNPs_after_plink_ld.5.no-tra.snp-sv.missing-filter.txt -n 18624 -o ~/projects/sv_nams/data/subset-NAM-snps/SNPs-to-keep_random-missing-filter_chr5.txt
shuf SNPs_after_plink_ld.6.no-tra.snp-sv.missing-filter.txt -n 16421 -o ~/projects/sv_nams/data/subset-NAM-snps/SNPs-to-keep_random-missing-filter_chr6.txt
shuf SNPs_after_plink_ld.7.no-tra.snp-sv.missing-filter.txt -n 13169 -o ~/projects/sv_nams/data/subset-NAM-snps/SNPs-to-keep_random-missing-filter_chr7.txt
shuf SNPs_after_plink_ld.8.no-tra.snp-sv.missing-filter.txt -n 16273 -o ~/projects/sv_nams/data/subset-NAM-snps/SNPs-to-keep_random-missing-filter_chr8.txt
shuf SNPs_after_plink_ld.9.no-tra.snp-sv.missing-filter.txt -n 14818 -o ~/projects/sv_nams/data/subset-NAM-snps/SNPs-to-keep_random-missing-filter_chr9.txt
shuf SNPs_after_plink_ld.10.no-tra.snp-sv.missing-filter.txt -n 14424 -o ~/projects/sv_nams/data/subset-NAM-snps/SNPs-to-keep_random-missing-filter_chr10.txt


# subsample random SNPs based on SNPs with very high or with very low LD
for chr in {1..10}; do
  qsub -v CHR=$chr ~/projects/sv_nams/scripts/subsample_high-low_ld.sh
done

# # results: 152,408 snps in high ld
# 24230 SNPs-to-keep_subsample-high-ld_chr1.missing-filter.txt
# 11885 SNPs-to-keep_subsample-high-ld_chr2.missing-filter.txt
# 16699 SNPs-to-keep_subsample-high-ld_chr3.missing-filter.txt
# 14540 SNPs-to-keep_subsample-high-ld_chr4.missing-filter.txt
# 16875 SNPs-to-keep_subsample-high-ld_chr5.missing-filter.txt
# 14793 SNPs-to-keep_subsample-high-ld_chr6.missing-filter.txt
# 12110 SNPs-to-keep_subsample-high-ld_chr7.missing-filter.txt
# 14928 SNPs-to-keep_subsample-high-ld_chr8.missing-filter.txt
# 13450 SNPs-to-keep_subsample-high-ld_chr9.missing-filter.txt
# 12898 SNPs-to-keep_subsample-high-ld_chr10.missing-filter.txt

# # results: 169,793 snps in low ld
# 27830 SNPs-to-keep_subsample-low-ld_chr1.missing-filter.txt
# 13358 SNPs-to-keep_subsample-low-ld_chr2.missing-filter.txt
# 18647 SNPs-to-keep_subsample-low-ld_chr3.missing-filter.txt
# 16229 SNPs-to-keep_subsample-low-ld_chr4.missing-filter.txt
# 18624 SNPs-to-keep_subsample-low-ld_chr5.missing-filter.txt
# 16421 SNPs-to-keep_subsample-low-ld_chr6.missing-filter.txt
# 13169 SNPs-to-keep_subsample-low-ld_chr7.missing-filter.txt
# 16273 SNPs-to-keep_subsample-low-ld_chr8.missing-filter.txt
# 14818 SNPs-to-keep_subsample-low-ld_chr9.missing-filter.txt
# 14424 SNPs-to-keep_subsample-low-ld_chr10.missing-filter.txt
```

> The reason why the subset with SNPs in high LD with an SV had less SNPs than 169k is that plink doesn't compute LD for markers that have minimum allele frequency less than 0.05 or that are monomorphic. Since the difference is minimal, it is unlikely to interfere with the GWAS.

I ploted and summarized the R2 distribution for each subset.

```bash
cd ~/projects/sv_nams/analysis/reseq_snps_projection2

for chr in {1..10}; do
  qsub -v CHR=$chr ~/projects/sv_nams/scripts/distribution_snp-sv_ld.sh
done
```

Now, I have to filter the hapmap files to include only the SNPs from each subset, and then merged each hapmap dataset separately.

```bash
cd ~/projects/sv_nams/analysis/reseq_snps_projection2

for chr in {1..10}; do
  qsub -v CHR=$chr ~/projects/sv_nams/scripts/subset_hmp_based_on_ld.sh
done

# quick check results
for subset in snps-low-ld-sv snps-high-ld-sv snps-random; do
  wc -l NAM_rils_subset_SNPs.chr-*.$subset.hmp.txt
  echo " "
done

# merge all chromosomes
cd ~/projects/sv_nams/analysis/reseq_snps_projection2

# low ld
cp NAM_rils_subset_SNPs.chr-1.snps-low-ld-sv.hmp.txt NAM_rils_subset_SNPs.snps-low-ld-sv.hmp.txt
for chr in {2..10}; do
  echo $chr
  sed 1d NAM_rils_subset_SNPs.chr-$chr.snps-low-ld-sv.hmp.txt >> NAM_rils_subset_SNPs.snps-low-ld-sv.hmp.txt
done

# high ld
cp NAM_rils_subset_SNPs.chr-1.snps-high-ld-sv.hmp.txt NAM_rils_subset_SNPs.snps-high-ld-sv.hmp.txt
for chr in {2..10}; do
  echo $chr
  sed 1d NAM_rils_subset_SNPs.chr-$chr.snps-high-ld-sv.hmp.txt >> NAM_rils_subset_SNPs.snps-high-ld-sv.hmp.txt
done

# random snps
cp NAM_rils_subset_SNPs.chr-1.snps-random.hmp.txt NAM_rils_subset_SNPs.snps-random.hmp.txt
for chr in {2..10}; do
  echo $chr
  sed 1d NAM_rils_subset_SNPs.chr-$chr.snps-random.hmp.txt >> NAM_rils_subset_SNPs.snps-random.hmp.txt
done
```

I also created a hapmap only with SV subset to be used in GWAS comparisons.

```bash
cd ~/projects/sv_nams/analysis/reseq_snps_projection2

for chr in {1..10}; do
  qsub -v CHR=$chr ~/projects/sv_nams/scripts/subset_hmp_based_on_ld_svs.sh
done

# merge all chromosomes
cp NAM_rils_subset_SVs.chr-1.hmp.txt NAM_rils_subset_SVs.hmp.txt
for chr in {2..10}; do
  echo $chr
  sed 1d NAM_rils_subset_SVs.chr-$chr.hmp.txt >> NAM_rils_subset_SVs.hmp.txt
done
```

After creating all subsets, I performed some QC to make sure the data is ready for GWAS. I ploted the distribution of SNPs along the chromosome to see where the markers are in the genome, ploted the distribution of missing data for each dataset, and also the distribution of R2 values between SNPs and SVs.

```bash
cd ~/projects/sv_nams/analysis/reseq_snps_projection2

# distribution SNPs (or SNPs + SV) along chromosomes
Rscript ~/projects/sv_nams/scripts/distribution_snps-svs_chrom.R ld/subset_high-ld-snps/SNPs-kept_chr1.txt \
                                                                 ld/subset_high-ld-snps/distribution_snps_chrom_high.png

Rscript ~/projects/sv_nams/scripts/distribution_snps-svs_chrom.R ld/subset_low-ld-snps/SNPs-kept_chr1.txt \
                                                                 ld/subset_low-ld-snps/distribution_snps_chrom_low.png

Rscript ~/projects/sv_nams/scripts/distribution_snps-svs_chrom.R ld/subset_random-snps/SNPs-kept_chr1.txt \
                                                                 ld/subset_random-snps/distribution_snps_chrom_random.png


# histogram percent missing data for each subset
run_pipeline.pl -Xmx40g -importGuess NAM_rils_subset_SNPs.snps-high-ld-sv.hmp.txt \
                -GenotypeSummaryPlugin -endPlugin \
                -export ld/subset_high-ld-snps/tassel_summary

run_pipeline.pl -Xmx40g -importGuess NAM_rils_subset_SNPs.snps-low-ld-sv.hmp.txt \
                -GenotypeSummaryPlugin -endPlugin \
                -export ld/subset_low-ld-snps/tassel_summary

run_pipeline.pl -Xmx40g -importGuess NAM_rils_subset_SNPs.snps-random.hmp.txt \
                -GenotypeSummaryPlugin -endPlugin \
                -export ld/subset_random-snps/tassel_summary

run_pipeline.pl -Xmx40g -importGuess NAM_rils_subset_SVs.hmp.txt \
                -GenotypeSummaryPlugin -endPlugin \
                -export ld/tassel_summary_sv

Rscript ~/projects/sv_nams/scripts/qc_tassel_summary.R ld/subset_high-ld-snps/tassel_summary3.txt \
                                                       ld/subset_high-ld-snps/missing_snps_high.png

Rscript ~/projects/sv_nams/scripts/qc_tassel_summary.R ld/subset_low-ld-snps/tassel_summary3.txt \
                                                       ld/subset_low-ld-snps/missing_snps_low.png

Rscript ~/projects/sv_nams/scripts/qc_tassel_summary.R ld/subset_random-snps/tassel_summary3.txt \
                                                       ld/subset_random-snps/missing_snps_random.png

Rscript ~/projects/sv_nams/scripts/qc_tassel_summary.R ld/tassel_summary_sv3.txt \
                                                       ld/missing_svs_subset.png

# plot distribution R2
Rscript ~/projects/sv_nams/scripts/distribution_snps-LD-svs_all-chr.R ld/subset_high-ld-snps \
                                                                      ld/subset_high-ld-snps/dist-LD_SNPs-SVs_high.png

Rscript ~/projects/sv_nams/scripts/distribution_snps-LD-svs_all-chr.R ld/subset_low-ld-snps \
                                                                      ld/subset_low-ld-snps/dist-LD_SNPs-SVs_low.png

Rscript ~/projects/sv_nams/scripts/distribution_snps-LD-svs_all-chr.R ld/subset_random-snps \
                                                                      ld/subset_random-snps/dist-LD_SNPs-SVs_random.png
```

Finally, I uploaded the 4 subsets (3 SNP subsets and 1 SV subset) to Cyverse:

```bash
# go to data folder of the project
cd ~/projects/sv_nams/analysis/reseq_snps_projection2

# compress for faster upload
gzip NAM_rils_subset_SNPs.snps-high-ld-sv.hmp.txt
gzip NAM_rils_subset_SNPs.snps-low-ld-sv.hmp.txt
gzip NAM_rils_subset_SNPs.snps-random.hmp.txt
gzip NAM_rils_subset_SVs.hmp.txt

# rename files for consistency with versions of projected snps and svs
cp NAM_rils_subset_SNPs.snps-high-ld-sv.hmp.txt.gz NAM_rils_subset_SNPs.snps-high-ld-sv.v7.hmp.txt.gz
cp NAM_rils_subset_SNPs.snps-low-ld-sv.hmp.txt.gz NAM_rils_subset_SNPs.snps-low-ld-sv.v7.hmp.txt.gz
cp NAM_rils_subset_SNPs.snps-random.hmp.txt.gz NAM_rils_subset_SNPs.snps-random.v7.hmp.txt.gz
cp NAM_rils_subset_SVs.hmp.txt.gz NAM_rils_subset_SVs.v7.hmp.txt.gz

# log in to cyverse
iinit
# go to cyverse shared folder to download data
icd /iplant/home/shared/NAM/Misc
# check if files match what Arun described
ils
# upload data
iput -K NAM_rils_subset_SNPs.snps-high-ld-sv.v7.hmp.txt.gz
iput -K NAM_rils_subset_SNPs.snps-low-ld-sv.v7.hmp.txt.gz
iput -K NAM_rils_subset_SNPs.snps-random.v7.hmp.txt.gz
iput -K NAM_rils_subset_SVs.v7.hmp.txt.gz
# exit iRods
iexit full
```



# Full SNP and SV dataset

Dr. Tingting requested datasets containing information about all SNPs and all SVs separately (and split into 10 chromosomes each). I already have the files ready for SVs, but for SNP files I will need to quickly filter the final file with both SNPs and SVs. Then I just need compress the files and upload them into Cyverse.

```bash
cd ~/projects/sv_nams/analysis/reseq_snps_projection2

# svs
for chr in {1..10}; do
  echo $chr
  gzip -c NAM_rils_projected-SVs-only.all-RILs.chr-$chr.hmp.txt > NAM_rils_projected-SVs-only.all-RILs.chr-$chr.v7.hmp.txt.gz
done

# snps
for chr in {1..10}; do
  qsub -v CHR=$chr ~/projects/sv_nams/scripts/split_reseq-snps_by_chr.sh
done

# log in to cyverse
iinit
# go to cyverse shared folder to download data
icd /iplant/home/shared/NAM/Misc
# check if files match what Arun described
ils
# upload SV dataset
for chr in {1..10}; do
  iput -K NAM_rils_projected-SVs-only.all-RILs.chr-$chr.v7.hmp.txt.gz
done
# upload SNP datasets
for chr in {1..10}; do
  iput -K NAM_rils_projected-reseq-SNPs-only.all-RILs.chr-$chr.v7.hmp.txt.gz
done
# exit iRods
iexit full
```

All files were moved to Cyverse folder: `/iplant/home/shared/NAM/Misc/NAM-SV-projected-V7`
