#!/bin/bash
nam=$1
index=$2
vcf=$3
fasta=$4
base="/work/LAS/mhufford-lab/arnstrm/Canu_1.8/genetic_maps"
om1_splitOnemap.sh results_onemap.txt
for f in *.csv; do
   om2_CheckSplitsConsistency.sh $f;
done
for f in *.csv; do
   om3_FixLinkageGroupName.sh $f;
done
om4_MakeMarkerFile.sh ${nam}
om5_ExtractAndMapMarkers.sh ${index} ${nam}_onemap-cM.csv ${vcf} ${nam}
#singularity shell --bind $PWD ${base}/jcvi.simg ${base}/99_runALLMAPS.sh ${nam}_mapped_onemap-cM.csv ${fasta}
