#!/bin/bash
source /work/LAS/mhufford-lab/arnstrm/miniconda3/etc/profile.d/conda.sh
conda activate bionano
#source activate bionano
ml r-devtools/1.12.0-py2-r3.4-47jzlan
ml rm perl
ml rm python
Rscript -e '.libPaths("~/R/x86_64-pc-linux-gnu-library/3.4/")'
EBROOTBIONANOSOLVE="/work/LAS/mhufford-lab/arnstrm/bionano/Solve3.4_06042019a"
rcmap="/work/LAS/mhufford-lab/arnstrm/stiff_stalks/Zea_mays_var_B73_V4_release/enz2/zm_B73_v4_all_chrs_un_ctgs-chr-only_GCTCTTC_0kb_0labels.cmap"
qcmap=$1
python $EBROOTBIONANOSOLVE/Pipeline/06042019/runCharacterize.py \
-t $EBROOTBIONANOSOLVE/RefAligner/8949.9232rel/RefAligner \
-r ${rcmap} \
-q ${qcmap} \
-p $EBROOTBIONANOSOLVE/Pipeline/06042019 \
-a $EBROOTBIONANOSOLVE/RefAligner/8949.9232rel/optArguments_nonhaplotype_noES_DLE1_saphyr.xml \
-n 32 > B73v4-${qcmap%.*}.stats

python $EBROOTBIONANOSOLVE/Pipeline/06042019/runSV.py \
-t $EBROOTBIONANOSOLVE/RefAligner/8949.9232rel/RefAligner \
-r ${rcmap} \
-q ${qcmap} \
-o B73v4-${qcmap%.*} \
-E alignref/${qcmap%.*}.errbin \
-a $EBROOTBIONANOSOLVE/RefAligner/8949.9232rel/optArguments_nonhaplotype_noES_DLE1_saphyr.xml \
-j 32 \
-T 32 \
