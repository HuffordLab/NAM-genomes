#!/bin/bash
# run this on right files
# merges them to one file
if [ "$#" -ne 1 ] ; then
echo "please provide a NAM name for the marker file";
exit 0;
fi
base=$1
rm ${base}_onemap-cM.csv 2>/dev/null;
cat *_right.csv | grep -v "^chr" >> ${base}_onemap-cM.csv;
