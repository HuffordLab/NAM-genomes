#!/bin/bash
file=$1
echo -e "PS\tWorkflow\tGenes"
awk '{print $1"\tWF1\t"$2}' $file |grep -v "^PS"
awk '{print $1"\tWF1-Filtered\t"$3}' $file | grep -v "^PS"
awk '{print $1"\tWF2\t"$4}' $file |grep -v "^PS"
awk '{print $1"\tWF2-Filtered\t"$5}' $file |grep -v "^PS"
awk '{print $1"\tTotal\t"$6}' $file |grep -v "^PS"
