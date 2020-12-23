#!/bin/bash
module load picard
module load gatk
file="$1"
module load r-geneplotter/1.58.0-py2-r3.5-ghvy7k7
ref="/work/LAS/mhufford-lab/arnstrm/NAM/rnaseq/00_B73index/B73.fasta"
#TZI-8_8DAS_root_MN27011.bam
id=$(echo ${file%.*} |cut -f 1,4 -d "_")
mn=$(echo $id | rev |cut -c 1-7 |rev)
name=$(echo ${file%.*} |cut -f 2,3 -d "_")

rg="${id}"
machine="NovaSeq"
platform="illumina"
sample="${id}"
lib="${mn}"

echo $rg $lib $platform $machine  $sample
picard BuildBamIndex INPUT=$file
picard AddOrReplaceReadGroups \
      I=$file \
      O=${file%.*}_rg.bam \
      SO=coordinate \
      RGID=${rg} \
      RGLB=${lib} \
      RGPL=${platform} \
      RGPU=${machine} \
      RGSM=$sample

picard MarkDuplicates \
     INPUT=${file%.*}_rg.bam \
     OUTPUT=${file%.*}_rg-md.bam \
     METRICS_FILE=${file%.*}_markduplicates_metrics.txt \
     VALIDATION_STRINGENCY=SILENT \
     OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
     CREATE_INDEX=true \
     TMP_DIR=$TMPDIR

gatk SplitNCigarReads \
   --reference ${ref} \
   --create-output-bam-index true \
   --input ${file%.*}_rg-md.bam \
   --output ${file%.*}_split.bam
