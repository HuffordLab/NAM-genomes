#!/bin/bash
if [ "$#" -ne 4 ] ; then
echo "please provide:"
echo -e "\t\t(1) path for hisat2 indexed NAM genome (without ht2 extension)"
echo -e "\t\t(2) csv marker file checked, processed and merged"
echo -e "\t\t(3) VCF file to obtain co-ordinates for the markers"
echo -e "\t\t(4) NAM line name for naming the output file"
echo "";
echo "./04_onemap_ExtractAndMapMarkers.sh <index_path> <csv_file> <vcf_file> <NAM-name>" ;
echo "";
exit 0;
fi
module load bedtools2
module load hisat2
module load samtools
module load bamtools
module load star
B73v4="/work/LAS/mhufford-lab/arnstrm/Canu_1.8/required-files/B73.fasta"
# full path for ht2 files
index="$1"
# full path for markers csv file
csv="$2"
# NAM genome
#genome="$3"
# vcf file for locations
vcf="$3"
# full nam name eg Ki3
full="$4"
# convert the genetic map to tab and only keep SNP name and genetic distance
sed 's/,/\t/g' $csv | grep -v "chr" | cut -f 2,3 > marker_distnace.txt
# get the B73.V4 locations for SNPs from the subsetted VCF file
grep -Fw -f <(cut -f 1 marker_distnace.txt) ${vcf} | awk '{print $3"\t"$1"\t"$2}' > marker_locations.txt
# merge locations with distance and create a bed file for extracting marker sequence
# uses 50bp up and down from the SNP location
awk 'BEGIN{OFS=FS="\t"}FNR==NR{a[$1]=$2 FS $3;next}{ print $0, a[$1]}' marker_distnace.txt marker_locations.txt |\
  sed 's/ //g' |\
  awk 'BEGIN{OFS=FS="\t"}{print $2,$3-50,$3+50,$1,$4,"+"}' > markers_locations_distance.bed
# use bedtools to extract marker sequence
bedtools getfasta -fi ${B73v4} -fo ${full}_markers.fasta -bed markers_locations_distance.bed -name
# save marker sequnce as variable
markers=${full}_markers.fasta
# generate index for the genome
#hisat2-build ${genome} ${base}_index
# map marker sequence to the indexed genome
#hisat2 -p 12 --mp 1,1 --no-softclip -f -x ${index} -U ${markers}  1> ${nam}_mapped.sam 2> mapping_stats.txt
#STAR --runMode alignReads --runThreadN 16 --genomeDir ${index} --outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 --outFileNamePrefix ${markers%.*}-star --readFilesIn ${markers}
# sort and covert to bam
#hisat2 -p 12 --mp 1,1 --no-softclip -f -x $index -U ${markers}  -S ${markers%.*}_mapped.sam &> mapping_stats.txt
# convert sam to bam and then sort
#samtools view -b -o ${markers%.*}_mapped.bam ${markers%.*}-starAligned.out.sam
#samtools sort -o ${markers%.*}_mapped_sorted.bam ${markers%.*}_mapped.bam
# save the marker name and location they map
#bedtools bamtobed -i ${markers%.*}_mapped_sorted.bam |  awk '{print $4"\t"$1"\t"$2}' > ${markers%.*}_part1.txt
# get the linkage group from the genetic map
hisat2 -p 12 --mp 1,1  \
    --no-spliced-alignment \
    --end-to-end \
    --rdg 10000,10000 \
    --rfg 10000,10000 \
    -f -x ${index} -U ${markers}  1> ${nam}_mapped.sam 2> mapping_stats.txt
samtools view -b -o ${nam}_mapped.bam ${nam}_mapped.sam
samtools sort -o ${nam}_sorted.bam ${nam}_mapped.bam
samtools view -q 30 -b ${nam}_sorted.bam > ${nam}_sorted-q30.bam
bamtools filter -tag XM:0 -in ${nam}_sorted-q30.bam -out ${nam}_sorted_noMismatch.bam
bedtools bamtobed -i ${nam}_sorted_noMismatch.bam |  awk '{print $4"\t"$1"\t"$2}' > ${nam}_part1.txt


awk 'BEGIN{OFS=FS=","} {print $2,$1,$3}' $csv | sed 's/,/\t/g' > ${markers%.*}_part2.txt
# write a csv file merging the part1 and part2 file
echo "Scaffold ID,scaffold position,LG,genetic position" > ${full}_mapped_onemap-cM.csv
awk 'BEGIN{OFS=FS="\t"}FNR==NR{a[$1]=$2 FS $3;next}{ print $0, a[$1]}' ${markers%.*}_part2.txt ${markers%.*}_part1.txt | \
  sed 's/ //g' | \
  cut -f 2- | \
  sed 's/\t/,/g' >> ${full}_mapped_onemap-cM.csv
