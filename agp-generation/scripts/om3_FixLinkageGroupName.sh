#!/bin/bash
# correct the linkage group to the actual chromosome name based on B73 naming
if [ "$#" -ne 1 ] ; then
echo "please provide a NAM name for the marker file";
exit 0;
fi
file="$1"
correct=$(cut -f 2 -d "," $file | cut -f 1 -d "_" | sed 's/S//g' | sort | uniq -c | awk '{print $2"\t"$1}' | sort -k 2 -n | tail -n 1 |cut -f 1);
current=$(cut -f 1 -d "," $file | sort | uniq )
awk -v x=${correct} 'BEGIN{OFS=FS=","} $1=x'  $file > ${file%.*}_right.csv
echo "the linkage group was listed as $current but was changed to $correct"


