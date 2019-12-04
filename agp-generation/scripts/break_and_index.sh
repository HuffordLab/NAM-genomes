#!/bin/bash
module use /work/GIF/software/modules
module load GIF2/UCSC
NAM=$(basename $(dirname $(pwd)))
file=$1
faSplit size $file 100000000 ${file%.*}-broken -oneFile
module load hisat2
hisat2-build ${file%.*}-broken.fa $NAM
