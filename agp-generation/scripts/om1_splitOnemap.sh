#!/bin/bash
# run this onemap result file
# splits them based on linkage groups to check consistency
if [ "$#" -ne 1 ] ; then
echo "please provide the onemap results file";
fi

file="$1"
base=$(basename ${file%.*} |cut -f 1 -d "_")
min=$(cut -f 1 -d " " ${file} | awk 'NF>0' | sort -n | uniq | head -n 1)
max=$(cut -f 1 -d " " ${file} | awk 'NF>0' | sort -n | uniq | tail -n 1)
for i in $(seq ${min} ${max}); do
awk -v x=$i '$1==x' $file | sed 's/ /,/g' > ${base}_onemap_chr${i}.csv;
done
