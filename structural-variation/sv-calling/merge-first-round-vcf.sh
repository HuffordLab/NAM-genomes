#!/bin/bash
ls *.vcf > vcf.fofn
SURVIVOR merge vcf.fofn 1000 1 1 -1 -1 -1 survivor-merged_1kbpdist_r1.vcf
