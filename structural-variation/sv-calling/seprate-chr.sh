#!/bin/bash
genome=$1
num=$(grep -n ">" $genome | head -n 11 |tail -n 1 |cut -f 1 -d ":" |awk '{print $1-1}' )
head -n $num $genome > ${genome%.*}-chr-only.fa

