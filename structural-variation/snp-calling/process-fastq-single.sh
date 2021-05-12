
#!/bin/bash
module load picard
module load bwa
module load samtools
module load parallel

PICARD_HOME="/work/LAS/mhufford-lab/arnstrm/newNAM/analyses/l-gatk-illumina"
REF="/work/LAS/mhufford-lab/arnstrm/newNAM/analyses/l-gatk-illumina/b-genome/B73.PLATINUM.pseudomolecules-v1.fasta"
ulimit -c unlimited
PLT="illumina"
TDATE=$(date '+%Y-%m-%d %H:%M:%S' |sed 's/ /T/g')

# check
echo $PLT
echo $TDATE

R1="${1}"
R2="${2}"

# convert fastq to sam and add readgroups
#parallel -j 4 --link "java -Djava.io.tmpdir=$TMPDIR -Xmx150G -jar $PICARD_HOME/picard.jar FastqToSam FASTQ={1} FASTQ2={2} OUTPUT={1.}_fastqtosam.bam READ_GROUP_NAME={1.} SAMPLE_NAME={1.}_name LIBRARY_NAME={1.}_lib PLATFORM_UNIT=${PLT} PLATFORM=illumina SEQUENCING_CENTER=ISU RUN_DATE=${TDATE}" ::: ${R1} ::: ${R2}
# mark adapters
#parallel -j 4 "java -Djava.io.tmpdir=$TMPDIR -Xmx150G -jar $PICARD_HOME/picard.jar MarkIlluminaAdapters I={.}_fastqtosam.bam O={.}_markilluminaadapters.bam M={.}_markilluminaadapters_metrics.txt TMP_DIR=${TMPDIR}" ::: ${R1}
# convert bam back to fastq for mapping
#parallel -j 4 "java -Djava.io.tmpdir=$TMPDIR -Xmx150G -jar $PICARD_HOME/picard.jar SamToFastq I={.}_markilluminaadapters.bam FASTQ={.}_samtofastq_interleaved.fq CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true TMP_DIR=${TMPDIR}" ::: ${R1}
# mapping reads to indexed genome
#parallel -j 1 "bwa mem -M -t 35 -p $REF {.}_samtofastq_interleaved.fq | samtools view -buS - > {.}_bwa_mem.bam" ::: ${R1}
# merging alignments
parallel -j 4 "java -Djava.io.tmpdir=$TMPDIR -Xmx150G -jar $PICARD_HOME/picard.jar MergeBamAlignment R=$REF UNMAPPED_BAM={.}_fastqtosam.bam ALIGNED_BAM={.}_bwa_mem.bam O={.}_snippet_mergebamalignment.bam CREATE_INDEX=true ADD_MATE_CIGAR=true CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS TMP_DIR=${TMPDIR}" ::: ${R1}
# mark duplicates
parallel -j 4 "java -Djava.io.tmpdir=$TMPDIR -Xmx150G -jar $PICARD_HOME/picard.jar MarkDuplicates INPUT={.}_snippet_mergebamalignment.bam OUTPUT={.}_final.bam METRICS_FILE={.}_mergebamalignment_markduplicates_metrics.txt OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 CREATE_INDEX=true TMP_DIR=$TMPDIR" ::: ${R1}
echo >&2 "ALL DONE!"
#parallel "rm {.}_fastqtosam.bam {.}_markilluminaadapters.bam {.}_samtofastq_interleaved.fq {.}_bwa_mem.bam {.}_snippet_mergebamalignment.bam {.}_snippet_mergebamalignment.bai" ::: ${R1}
