#!/bin/bash

dir=$1

cd $dir
rm ${dir}.cmds
for f in error_corrected_reads/*.gz; do
echo "./runNGMLR.sh $f";
done >> ${dir}.cmds
makeSLURMp.py 36 ${dir}.cmds;
for sub in *.sub; do
sed -i 's/parallel -j 1/parallel -j 9/g' ${sub};
done
