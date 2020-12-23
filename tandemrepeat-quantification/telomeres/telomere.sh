#!/bin/bash
#This script find the trf coordinates on the small and long arm of each chromosome of the input genome 
module load trf
for file in *fasta; do
trf $file 2 7 7 80 10 50 500 -f -d -m -h
awk '{if ($3==7) print}' $file.2.7.7.80.10.50.500.dat | head -2 | tail -1 | grep "CCCTAAA" | awk '{print $1"\t"$2}' > $file.small_arm_coordinates
awk '{if ($3==7) print}' $file.2.7.7.80.10.50.500.dat | grep "TTTAGGG" | tail -1 | awk '{print $1"\t"$2}' > $file.long_arm_coordinates
done
