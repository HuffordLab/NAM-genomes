#!/bin/bash
module load star
base=$(basename $(pwd))
genome="$1"
indexdir="$(basename ${genome%.*})-star"
mkdir -p $indexdir
STAR \
--runMode genomeGenerate \
--runThreadN 36 \
--genomeDir ${indexdir} \
--genomeFastaFiles ${genome}
