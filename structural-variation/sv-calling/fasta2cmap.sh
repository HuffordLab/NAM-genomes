#!/bin/bash
source /work/LAS/mhufford-lab/arnstrm/miniconda/etc/profile.d/conda.sh
conda activate bionano
ml r-devtools/1.12.0-py2-r3.4-47jzlan
ml rm perl
ml rm python
Rscript -e '.libPaths("~/R/x86_64-pc-linux-gnu-library/3.4/")'
EBROOTBIONANOSOLVE="/work/LAS/mhufford-lab/arnstrm/PanAnd/Solve3.4_06042019a"
genome=$1
perl $EBROOTBIONANOSOLVE/HybridScaffold/06042019/scripts/fa2cmap_multi_color.pl -i $genome -e CTTAAG 1 -e CACGAG 2
