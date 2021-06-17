#!/bin/bash
module load samtools
dir=$1
cd $dir
mkdir -p samfiles bamfiles sortedbam
for sam in *.sam; do
ibam=${sam%.*}.bam
sbam=${sam%.*}_sorted.bam
samtools view --threads 36 -b -o $ibam $sam
exit_status=$?
if [ $exit_status -eq 1 ]; then
    echo -e "$sam\t ERROR"
fi
samtools sort -o $sbam -T ${sam%.*}_temp --threads 36 $ibam
mv $sam ./samfiles/
mv $ibam ./bamfiles/
mv $sbam ./sortedbam/
done
cd ./sortedbam/
#module load bamtools
ls *.bam > input.fofn
samtools merge --threads 36 -O BAM -b input.fofn ../NAM_merged_ec_reads_mapped_to_B73.bam
#bamtools merge -list input.fofn -out ../NAM_merged_ec_reads_mapped_to_B73.bam
