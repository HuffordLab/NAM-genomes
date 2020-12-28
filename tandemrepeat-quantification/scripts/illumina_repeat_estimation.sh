#!/bin/bash

line="$1"

module load seqtk
module load BLAST/2.2.26-Linux_x86_64
module load BEDTools/2.28.0-foss-2018a

r1_fq=/scratch/jl03308/NAM_pancentromere/trimmed/Input/${line}/*_30X_R1_val_1.fq.gz
prefix=$(basename $r1_fq |cut -f1 -d ".")

mkdir -p /scratch/jl03308/NAM_pancentromere/analysis/repeat_cal/${line}
mkdir -p /scratch/jl03308/NAM_pancentromere/analysis/repeat_cal/${line}/illumina 
mkdir -p /scratch/jl03308/NAM_pancentromere/analysis/repeat_cal/${line}/genome 

cd /scratch/jl03308/NAM_pancentromere/analysis/repeat_cal/${line}/illumina
zcat $r1_fq | seqtk seq - -A > ${prefix}.fa
seqtk sample -s 11 ${prefix}.fa 0.2 > 3x_1.fa
formatdb -p F -i 3x_1.fa

blastall -p blastn -i 3x_1.fa -d /scratch/jl03308/pacbio_denovo_assembly/reference/knob180.fa -m 8 -b 5000 -F F > knob180.blast
blastall -p blastn -i 3x_1.fa -d /scratch/jl03308/pacbio_denovo_assembly/reference/TR-1.fa -m 8 -b 5000 -F F > TR-1.blast
blastall -p blastn -i 3x_1.fa -d /scratch/jl03308/pacbio_denovo_assembly/reference/CentC.fa -m 8 -b 5000 -F F > CentC.blast
blastall -p blastn -i 3x_1.fa -d /scratch/jl03308/pacbio_denovo_assembly/reference/rDNA_intergenic_spacer.fa -m 8 -b 5000 -F F > rDNA_intergenic_spacer.blast
blastall -p blastn -i 3x_1.fa -d /scratch/jl03308/pacbio_denovo_assembly/reference/subtelomeric_4-12-1.fa -m 8 -b 5000 -F F > subtelomeric_4-12-1.blast
blastall -p blastn -i 3x_1.fa -d /scratch/jl03308/NAM_Canu1.8_verified_version_1/annotation/pan_TE_v1.0/centromere_TEs/CRMs.fa -m 8 -b 5000 -F F > CRMs.blast

#coverage=2.86

fc=$(cat /scratch/jl03308/NAM_Canu1.8_verified_version_1/pseudomolecules/NAM.fc.genome.size |awk -v var1=$line '{if($1==var1){print$2}}')
rn=$(cat 3x_1.fa | grep -c ">")
coverage=$(echo $rn*150/$fc/1000000000 |bc -l)

for file in *.blast
do 
total=$(cat $file |awk '{if($4>=30){print$0}}'|cut -f1,7,8|bedtools sort -i - |bedtools merge -i - |awk '{print$0"\t"$3-$2}' |awk '{ SUM += $4} END { print SUM}')
adjust=$(echo $total/$coverage |bc)
echo -e $(basename $file |cut -f1 -d ".")"\t"$adjust >> $line.illumina
done

cd /scratch/jl03308/NAM_pancentromere/analysis/repeat_cal/${line}/genome
module load BLAST/2.2.26-Linux_x86_64

formatdb -p F -i /scratch/jl03308/NAM_Canu1.8_verified_version_1/pseudomolecules/${line}*pseudomolecules-v*.fasta
blastall -p blastn -d /scratch/jl03308/NAM_Canu1.8_verified_version_1/pseudomolecules/${line}*pseudomolecules-v*.fasta -i /scratch/jl03308/pacbio_denovo_assembly/reference/knob180.fa -m 8 -b 5000 > knob180.blast
blastall -p blastn -d /scratch/jl03308/NAM_Canu1.8_verified_version_1/pseudomolecules/${line}*pseudomolecules-v*.fasta -i /scratch/jl03308/pacbio_denovo_assembly/reference/TR-1.fa -m 8 -b 5000 > TR-1.blast
blastall -p blastn -d /scratch/jl03308/NAM_Canu1.8_verified_version_1/pseudomolecules/${line}*pseudomolecules-v*.fasta -i /scratch/jl03308/pacbio_denovo_assembly/reference/CentC.fa -m 8 -b 5000 > CentC.blast
blastall -p blastn -d /scratch/jl03308/NAM_Canu1.8_verified_version_1/pseudomolecules/${line}*pseudomolecules-v*.fasta -i /scratch/jl03308/pacbio_denovo_assembly/reference/rDNA_intergenic_spacer.fa -m 8 -b 5000 > rDNA_intergenic_spacer.blast
blastall -p blastn -d /scratch/jl03308/NAM_Canu1.8_verified_version_1/pseudomolecules/${line}*pseudomolecules-v*.fasta -i /scratch/jl03308/pacbio_denovo_assembly/reference/subtelomeric_4-12-1.fa -m 8 -b 5000 > subtelomeric_4-12-1.blast
blastall -p blastn -d /scratch/jl03308/NAM_Canu1.8_verified_version_1/pseudomolecules/${line}.*pseudomolecules-v*.fasta -i /scratch/jl03308/NAM_Canu1.8_verified_version_1/annotation/pan_TE_v1.0/centromere_TEs/CRMs.fa -m 8 -b 5000 > CRMs.blast

for file in *.blast
do 
echo $(basename $file |cut -f1 -d ".")
forward=$(cat $file |awk '{if($4>=30){print$0}}'|cut -f2,9,10 |awk '{if($3>$2){print$0}}' |bedtools sort -i - |bedtools merge -i - |awk '{print$0"\t"$3-$2}' |bedtools groupby -i - -g 1 -o sum -c 4 |cut -f2 |awk '{ SUM += $1} END { print SUM}')
reverse=$(cat $file |awk '{if($4>=30){print$0}}'|cut -f2,9,10 |awk '{if($3<$2){print$1"\t"$3"\t"$2}}' |bedtools sort -i - |bedtools merge -i - |awk '{print$0"\t"$3-$2}' |bedtools groupby -i - -g 1 -o sum -c 4 |cut -f2 |awk '{ SUM += $1} END { print SUM}')
total=$(echo $forward+$reverse |bc)
forward2=$(cat $file |awk '{if($4>=30&&$2~"chr"){print$0}}'|cut -f2,9,10 |awk '{if($3>$2){print$0}}' |bedtools sort -i - |bedtools merge -i - |awk '{print$0"\t"$3-$2}' |bedtools groupby -i - -g 1 -o sum -c 4 |cut -f2 |awk '{ SUM += $1} END { print SUM}')
reverse2=$(cat $file |awk '{if($4>=30&&$2~"chr"){print$0}}'|cut -f2,9,10 |awk '{if($3<$2){print$1"\t"$3"\t"$2}}' |bedtools sort -i - |bedtools merge -i - |awk '{print$0"\t"$3-$2}' |bedtools groupby -i - -g 1 -o sum -c 4 |cut -f2 |awk '{ SUM += $1} END { print SUM}')
total2=$(echo $forward2+$reverse2 |bc)
echo -e $(basename $file |cut -f1 -d ".")"\t"$total >> $line.genome.sum
echo -e $(basename $file |cut -f1 -d ".")"\t"$total2 >> $line.pseudo.sum
done

#paste $line.pseudo.sum $line.genome.sum ../illumina/$line.illumina.sum |cut -f1,2,4,6 |awk -v var1=$line '{print var1"\t"$0}' > /scratch/jl03308/NAM_pancentromere/analysis/repeat_cal/$line/$line.repeat.sum

paste <(sort -k1 $line.pseudo.sum) <(sort -k1 $line.genome.sum) <(sort -k1 ../illumina/$line.illumina) |cut -f1,2,4,6 |awk -v var1=$line '{print var1"\t"$0}' > /scratch/jl03308/NAM_pancentromere/analysis/repeat_cal/$line/$line.repeat.sum
