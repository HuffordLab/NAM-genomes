#!/bin/bash
if [ "$#" -ne 3 ] ; then
echo "please provide:"
echo -e "\t\t(1) AGP file"
echo -e "\t\t(2) NAM name"
echo -e "\t\t(3) Chromosome fasta file"
echo "";
echo "./02_gg_MapMarkers.sh <agp> <NAM-name> <pseudomolecules-fasta>" ;
echo "";
exit 0;
fi
agp=$1
nam=$2
chr=$3
# fix AGP
echo -ne "Writing new AGP file and versioning...\t"
grep "^#" $agp > ${nam}-pg_and_gg-based_v1.agp
grep "^chr" $agp | awk 'BEGIN{OFS=FS="\t"}$5=="W"{$6="scaf_"$6}{print}' >> ${nam}-pg_and_gg-based_v1.agp
grep -v "^chr" $agp |grep -v "^#" | awk 'BEGIN{OFS=FS="\t"}$5=="W"{$1="scaf_"$1}{$6="scaf_"$6}{print}' >> ${nam}-pg_and_gg-based_v1.agp
echo "DONE!"
# fix Scaf names
echo -ne "changing scaf names and versioning...\t"
module load bioawk
bioawk -c fastx '{print ">scaf_"$name"\n"$seq}' ${nam}.scaffolds.fasta | fold > ${nam}.scaffolds-v1.fasta
echo "DONE!"
# move scaf file
echo -ne "moving original scaf to star index folder...\t"
#mv ${nam}.scaffolds.fasta  /work/LAS/mhufford-lab/arnstrm/Canu_1.8/genetic_maps/${nam}/${nam}_star/
echo "DONE!"
# fix chr
echo -ne "changing names in chromosome file...\t"
sed '/>chr/!s/>/>scaf_/g' $chr > ${nam}.pseudomolecules-v1.fasta
echo "DONE!"
echo -ne "cleaning up...\t"
mv $chr $agp *.pdf ../
echo "DONE!"

