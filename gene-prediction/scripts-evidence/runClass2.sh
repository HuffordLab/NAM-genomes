#!/bin/bash
ml purge
PATH=$PATH:/work/LAS/mhufford-lab/arnstrm/programs/CLASS-2.1.7
ml gcc/7.3.0-xegsmw4
ml samtools
#module use /work/GIF/software/modules
#module load GIF/class2/2.1.7
bam="$1"
#SJ="$2"
#SJ="/ptmp/LAS/arnstrm/TransAssembly/SJ.all"
out=$(basename ${bam%.*})
/work/LAS/mhufford-lab/arnstrm/programs/CLASS-2.1.7/run_class.pl \
   -a $bam \
   -o ${out}_class.gtf \
   -p 36 \
   --verbose \
   --clean

