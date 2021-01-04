#!/bin/bash
bam_chip="$1"
bam_input="$2"

prefix=$(basename $bam_chip |cut -f1 -d ".")
binsize=5000
wdir=$(dirname $bam_chip)

cd $wdir
module load deepTools/3.3.1-intel-2019b-Python-3.7.4
#bamCoverage -b $bam -p 12 -bs $binsize -o ${prefix}.bedgraph -of bedgraph --minMappingQuality 20 --ignoreDuplicates
bamCompare -b1 $bam_chip -b2 $bam_input -p 12 -bs $binsize -o ${prefix}.ChIP_Input.RPKM.bedgraph -of bedgraph --minMappingQuality 20 --ignoreDuplicates --normalizeUsing RPKM --scaleFactorsMethod None
