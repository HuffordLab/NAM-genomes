#!/bin/bash
R1=$1
R2=$2
R1_base=$(echo $1 | sed 's/\.fastq//g')
R2_base=$(echo $2 | sed 's/\.fastq//g')
R1_name=${R1_base##*/}
R2_name=${R2_base##*/}
sub_base=${R1_name%??}

#10min for trimmomatic
cd /work/LAS/mhufford-lab/davehuf/NAMfoundersPrep
bash ./trim_pe.sh ${R1} ${R2}

#3.5 hours for bwa
bwa aln -t 64 GenomeData/ZmB73_RefGen_v2.allchr.mod.fasta ${R1_name}_paired.fq > ${R1_name}_paired.sai
bwa aln -t 64 GenomeData/ZmB73_RefGen_v2.allchr.mod.fasta ${R2_name}_paired.fq > ${R2_name}_paired.sai
bwa sampe GenomeData/ZmB73_RefGen_v2.allchr.mod.fasta ${R1_name}_paired.sai ${R2_name}_paired.sai ${R1_name}_paired.fq ${R2_name}_paired.fq > ${sub_base}_paired.sam

bwa aln -t 64 GenomeData/ZmB73_RefGen_v2.allchr.mod.fasta ${R1_name}_unpaired.fq > ${R1_name}_unpaired.sai
bwa samse GenomeData/ZmB73_RefGen_v2.allchr.mod.fasta ${R1_name}_unpaired.sai  ${R1_name}_unpaired.fq > ${sub_base}_1_unpaired.sam
bwa aln -t 64 GenomeData/ZmB73_RefGen_v2.allchr.mod.fasta ${R2_name}_unpaired.fq > ${R2_name}_unpaired.sai
bwa samse GenomeData/ZmB73_RefGen_v2.allchr.mod.fasta ${R2_name}_unpaired.sai ${R2_name}_unpaired.fq > ${sub_base}_2_unpaired.sam


samtools view -b -o ${sub_base}_paired.bam ${sub_base}_paired.sam
samtools view -b -o ${sub_base}_1_unpaired.bam ${sub_base}_1_unpaired.sam
samtools view -b -o ${sub_base}_2_unpaired.bam ${sub_base}_2_unpaired.sam

#five to six hours for add read groups for each bam
picard AddOrReplaceReadGroups INPUT=${sub_base}_paired.bam OUTPUT=${sub_base}_paired_RG.bam SORT_ORDER=coordinate RGID=paired RGLB=paired RGPL=illumina RGPU=${sub_base} RGSM=${sub_base} VALIDATION_STRINGENCY=LENIENT TMP_DIR=./
picard AddOrReplaceReadGroups INPUT=${sub_base}_1_unpaired.bam OUTPUT=${sub_base}_1_unpaired_RG.bam SORT_ORDER=coordinate RGID=unpaired_1 RGLB=unpaired_1 RGPL=illumina RGPU=${sub_base} RGSM=${sub_base} VALIDATION_STRINGENCY=LENIENT TMP_DIR=./
picard AddOrReplaceReadGroups INPUT=${sub_base}_2_unpaired.bam OUTPUT=${sub_base}_2_unpaired_RG.bam SORT_ORDER=coordinate RGID=unpaired_2 RGLB=unpaired_2 RGPL=illumina RGPU=${sub_base} RGSM=${sub_base} VALIDATION_STRINGENCY=LENIENT TMP_DIR=./

# three to four hours for merging bam files
picard MergeSamFiles INPUT=${sub_base}_paired_RG.bam INPUT=${sub_base}_1_unpaired_RG.bam INPUT=${sub_base}_2_unpaired_RG.bam OUTPUT=${sub_base}.bam SORT_ORDER=coordinate ASSUME_SORTED=false VALIDATION_STRINGENCY=LENIENT TMP_DIR=./
picard MarkDuplicates INPUT=${sub_base}.bam OUTPUT=${sub_base}.MD.bam METRICS_FILE=${sub_base}.MD.txt ASSUME_SORTED=true REMOVE_Duplicates=true VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=TRUE TMP_DIR=./
java -Xmx128g -jar /opt/rit/app/gatk/3.6/lib/GenomeAnalysisTK.jar -allowPotentiallyMisencodedQuals -I ${sub_base}.MD.bam -R GenomeData/ZmB73_RefGen_v2.allchr.mod.fasta -T RealignerTargetCreator -o ${sub_base}.forIndelRealigner.intervals
java -Xmx128g -jar /opt/rit/app/gatk/3.6/lib/GenomeAnalysisTK.jar -allowPotentiallyMisencodedQuals -I ${sub_base}.MD.bam -R GenomeData/ZmB73_RefGen_v2.allchr.mod.fasta -T IndelRealigner -targetIntervals ${sub_base}.forIndelRealigner.intervals -o ${sub_base}.IndelRealigned.bam
samtools index ${sub_base}.IndelRealigned.bam