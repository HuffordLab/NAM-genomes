#!/bin/bash
basedir="/work/LAS/mhufford-lab/arnstrm/Canu_1.8/genetic_maps"
nam=$1
index=$2
scaf=$3
for map in *.map; do
gg1_ProcessMaizeGDBmaps.sh $map $nam
done
gg2_MergeMarkers.sh $nam
gg3_MapMarkers.sh ${index} ${nam}_merged.fasta ${nam}_merged.txt ${nam}
singularity shell --bind $PWD ${basedir}/jcvi.simg ${basedir}/99_runALLMAPS.sh ${nam}_GG-mapped.csv ${scaf}
