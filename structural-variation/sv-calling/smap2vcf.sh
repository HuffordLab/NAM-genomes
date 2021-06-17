#!/bin/bash
source /work/LAS/mhufford-lab/arnstrm/miniconda3/etc/profile.d/conda.sh
conda activate bionano
ml r-devtools/1.12.0-py2-r3.4-47jzlan
ml rm perl
ml rm python
Rscript -e '.libPaths("~/R/x86_64-pc-linux-gnu-library/3.4/")'
EBROOTBIONANOSOLVE="/work/LAS/mhufford-lab/arnstrm/bionano/Solve3.4_06042019a"
nam=$1
#cd ${nam}/B73v5-*
pwd

python $EBROOTBIONANOSOLVE/VCFConverter/06042019/smap_to_vcf_v2.py \
-s ${nam}.smap \
-r ${nam}_r.cmap \
-x ${nam}.xmap \
-n ${nam} \
-o ${nam} \
-a "${nam}" \
-b False
