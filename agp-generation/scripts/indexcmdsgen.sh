#!/bin/bash
for f in */*_index/*.fasta; do
echo "00_build_index_hisat2.sh $(basename $f)" > $(dirname $f)/hisat.cmds;
cd $(dirname $f)
makeSLURMs.py 1 hisat.cmds;
cd /work/LAS/mhufford-lab/arnstrm/newNAM/analyses/d-genetic-maps
done

