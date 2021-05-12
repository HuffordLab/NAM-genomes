#!/bin/bash
source /work/LAS/mhufford-lab/arnstrm/miniconda/etc/profile.d/conda.sh
conda activate multiqc
#module load subread
GFF="/work/LAS/mhufford-lab/arnstrm/newNAM/analyses/f-rnaseq/index/GCA_000005005.6_B73_RefGen_v4_genomic.gtf"
base=$(basename $(pwd))
featureCounts \
    -a $GFF \
    -p \
    -o ${base}_counts.txt \
    -T 36 \
    -t exon \
    -g gene_id  *.bam

multiqc .
mv multiqc_report.html ${base}_qc.html
