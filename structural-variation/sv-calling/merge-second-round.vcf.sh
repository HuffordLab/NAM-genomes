#!/bin/bash
ls *_r2.vcf > vcf.fofn
SURVIVOR merge vcf.fofn 1000 -1 1 -1 -1 -1 survivor-merged_1kbpdist_r2.vcf
