#!/bin/bash

input="$1" #fq.gz file
read_num="$2"

ml seqtk/1.3-GCC-8.3.0
prefix=$(basename $input |cut -f1 -d "_")
wdir=$(dirname $input)
cd $wdir
seqtk sample -s100 $input $read_num | gzip > $wdir/${prefix}.subsample.fq.gz
