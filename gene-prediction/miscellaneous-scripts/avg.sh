#!/bin/bash
file=$1
awk '
   BEGIN { FS = "\t"; }
   NR == 1 { print; next; }
  { for (i = 3; i <= NF-1; i++) { sum[i] += $i; } }
  END {
    printf "\t";
    for (i = 3; i <= NF-1; i++) { printf "\t%.4g", sum[i] / (NR - 1); }
    print "";
  }
' $file
