#!/bin/bash
if [ "$#" -ne 2 ] ; then
echo "please provide:"
echo -e "\t\t(1) text file for markers dwonloaded from MaizeGDB, each CHR should be in separate file"
echo -e "\t\t(2) NAM line name for naming the output file"
echo "";
echo "./00_gg_ProcessMaizeGDBmaps.sh <maizegdb file> <NAM-name>" ;
echo "";
exit 0;
fi
# load needed modules
#module use /work/GIF/software/modules
module load cdbfasta
#module load GIF/maker
module load bioawk
module load bedtools2
# initialize variables
ref="/work/LAS/mhufford-lab/arnstrm/Canu_1.8/required-files/Zea_mays.AGPv3.22.dna.toplevel.fa"
# nam line name eg: Ki3
nam="$2"
# maizeGDB file (tsv format)
file="$1"
# read the file and choose the linkage group based on the markers mapping to B73 chr
chr=$(cut -f 5 $file |grep -v "C_v3" |sort |uniq -c |awk 'NF>1 {print $2"\t"$1}' |sort -k2,2 -rn |head -n 1 |cut -f 1)
# write a bed file for all markers with co-ordinates
cut -f 1,2,6 $file | awk 'NF==3' | sed 's/,//g' | awk -v x=$chr '{print x"\t"$3"\t"$3+100"\t"$1"\t"$2"\t+"}' > ${nam}_chr${chr}.bed
# extract the markers
bedtools getfasta -fi $ref -bed ${nam}_chr${chr}.bed -name > ${nam}_chr${chr}.fasta
# write a marker file with cM info
awk 'BEGIN{OFS=FS="\t"}{print $1"\t"$4"\t"$5}' ${nam}_chr${chr}.bed > ${nam}_chr${chr}.tsv
# markers without co-ordinates but with sequences instead
cut -f 1,2,6 $file | awk 'NF<3' | awk '{print $1}' > ${nam}_chr${chr}.notfound
# add them to the marker file with cM info and to the fasta file
grep -F -f ${nam}_chr${chr}.notfound $file | awk -v x=$chr 'NF==3 {print x"\t"$1"\t"$2}' >> ${nam}_chr${chr}.tsv
grep -F -f ${nam}_chr${chr}.notfound $file | awk 'NF==3 {print ">"$1"\n"$3}' >> ${nam}_chr${chr}.fasta
# another attempt to recover markers from Maggies data
cut -f 2 ${nam}_chr${chr}.tsv > ${nam}_chr${chr}.present
cut -f 1 $file | grep -v "Locus" | grep -Fwv -f ${nam}_chr${chr}.present > ${nam}_chr${chr}.absent
# Maggie's marker file
markers="/work/LAS/mhufford-lab/arnstrm/Canu_1.8/required-files/all_markers.txt"
# info extracted and written to a bed file
grep -Fwi -f ${nam}_chr${chr}.absent $markers | awk -v x=$chr '{print x"\t"$3"\t"$3+100"\t"$1"\t"$2"\t+"}' > ${nam}_chr${chr}.bed
# extract sequences
bedtools getfasta -fi $ref -bed ${nam}_chr${chr}.bed -name >> ${nam}_chr${chr}.fasta
awk 'BEGIN{OFS=FS="\t"}{print $1"\t"$4"\t"$5}' ${nam}_chr${chr}.bed >> ${nam}_chr${chr}.tsv
# remove intermediary files
rm ${nam}_chr${chr}.bed ${nam}_chr${chr}.notfound ${nam}_chr${chr}.present ${nam}_chr${chr}.absent
# print stats
ids=$(cut -f 2  ${nam}_chr${chr}.tsv | sort |uniq |wc -l)
seq=$(grep -c ">" ${nam}_chr${chr}.fasta)
tot=$(grep -v "Locus" $file |wc -l)
echo "For $nam Chr ${chr}, orignal file had ${tot} markers, but only $ids had information and $seq sequences were able to extract"
