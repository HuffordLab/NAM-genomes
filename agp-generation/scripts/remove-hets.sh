#!/bin/bash

het=$1
for markers in pangenome.csv goldengate.csv; do
while read line;  do
   awk -v x=$line -F"," '$1!=x' $markers  > temp
   mv temp $markers
done<$het
done
