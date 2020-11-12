#!/bin/bash
if [ "$#" -ne 4 ] ; then
echo "please provide:"
echo -e "\t\t(1) path for hisat2 indexed NAM genome (without ht2 extension)"
echo -e "\t\t(2) golder gate marker file, processed and merged"
echo -e "\t\t(3) information file for markers, with cM"
echo -e "\t\t(4) NAM line name for naming the output file"
echo "";
echo "./02_gg_MapMarkers.sh <index_path> <marker fasta> <info file txt> <NAM-name>" ;
echo "";
exit 0;
fi
# load modules
module load samtools
module load bamtools
module load bedtools2
module load hisat2
module load star
# variables
index="$1"
markers="$2"
nam="$4"
info="$3"
# map markers
hisat2 -p 12 --mp 1,1  \
    --no-spliced-alignment \
    --end-to-end \
    --rdg 10000,10000 \
    --rfg 10000,10000 \
    -f -x ${index} -U ${markers}  1> ${nam}_mapped.sam 2> mapping_stats.txt
#--no-softclip 
#STAR --runMode alignReads --runThreadN 16 --genomeDir ${index} --outFileNamePrefix ${nam}-star --readFilesIn ${markers}
#--outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 
# sort and covert to bam
samtools view -b -o ${nam}_mapped.bam ${nam}_mapped.sam
#${nam}-starAligned.out.sam
samtools sort -o ${nam}_sorted.bam ${nam}_mapped.bam
samtools view -q 30 -b ${nam}_sorted.bam > ${nam}_sorted-q30.bam
bamtools filter -tag XM:0 -in ${nam}_sorted-q30.bam -out ${nam}_sorted_noMismatch.bam

#bedtools bamtobed -i ${nam}_sorted.bam |  awk '{print $4"\t"$1"\t"$2}' > ${nam}_part1.txt
bedtools bamtobed -i ${nam}_sorted_noMismatch.bam |  awk '{print $4"\t"$1"\t"$2}' > ${nam}_part1.txt
# generate info
awk 'BEGIN{OFS=FS="\t"}{print $2,$1,$3}' ${info} > ${nam}_part2.txt
echo "Scaffold ID,scaffold position,LG,genetic position" > ${nam}_GG-mapped.csv
awk 'BEGIN{OFS=FS="\t"}FNR==NR{a[$1]=$2 FS $3;next}{ print $0, a[$1]}' ${nam}_part2.txt ${nam}_part1.txt | sed 's/ //g' | cut -f 2- |sed 's/\t/,/g' >> ${nam}_GG-mapped.csv
