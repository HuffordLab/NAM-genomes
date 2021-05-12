#!/bin/bash
input=$1
nam=${input%.pseudomolecules*}
new=${input#out_*_}
#awk -F '\t' -v nam="${nam}" '{ if ($3 == 100) { print $0 $nam } }' ${input} >> all100_${new}
head -n 1 ${input} | awk -F '\t' -v OFS="\t" -v n=${nam} '{print $0,n}' >> all_${new}

awk -F '\t' '{ if ( $3 == 100 && $9 == 100 ) {print $2,$7,$8,$1} }' ${input} >> QTL_markers_V5_coordinates.bed