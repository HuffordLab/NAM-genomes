#!/bin/bash
file=$1
cut -f 1,7-8 ${file} > ${file}.temp.bed
bedtools sort -i ${file}.temp.bed > ${file}.temp2.bed
bedtools merge -i ${file}.temp2.bed > ${file}.temp3.bed
cut -f 1 ${file}.temp3.bed |sort |uniq -c |awk '$1>1 {print $2"\t"$1}' > ${file}.duplicates.txt

