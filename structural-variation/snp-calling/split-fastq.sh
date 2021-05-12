#!/bin/bash
R1=$2
R2=$3
nam=$1
module load parallel
parallel <<FIL
zcat ${R1} | split -l 400000000 - ${nam}_R1_part_
zcat ${R2} | split -l 400000000 - ${nam}_R2_part_
FIL
mkdir -p gzdir
mv ${R1} ${R2} gzdir/
parallel "gzip {}" ::: ${nam}_R?_part_a?

