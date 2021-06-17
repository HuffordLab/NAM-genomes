#!/bin/bash

file=$1

bunzip2 $file
fasta-splitter.pl --n-parts 500 ${file%.*}
for f in  ${file%.*}-part*; do
gzip $f;
done
