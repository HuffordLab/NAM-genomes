#!/bin/bash

ChIP="$1"
Input="$2"
chrsize="$3"

prefix=$(basename $ChIP |cut -f1 -d ".")
wdir=$(dirname $ChIP)

cd $wdir

ml Miniconda3/4.7.10
source activate /home/jl03308/env-epic2
epic2 --control $Input --treatment $ChIP --mapq 20 --bin-size 5000 --gaps-allowed 0 --effective-genome-fraction 0.8 --chromsizes $chrsize --output ${prefix}.islands
conda deactivate
