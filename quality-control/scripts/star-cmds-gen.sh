#!/bin/bash
dir=$1
cd $dir
for r1 in *_R1.fq.gz; do
echo "../runSTARmap_round1.sh $r1";
done
echo "awk -f ../sjCollapseSamples.awk *_SJ.out.tab | sort -k1,1V -k2,2n -k3,3n > SJ.all";
for r1 in *_R1.fq.gz; do
echo "../runSTARmap_round2.sh $r1";
done
echo "../runSubreads.sh";
cd ..
