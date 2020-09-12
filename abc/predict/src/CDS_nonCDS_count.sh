snpcount=`zcat data/variants/B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.vcf.gz | grep -v "#" | grep -v scaf | wc -l`
echo "snp total $snpcount"

svcount=`cat data/variants/NAM_founders_SVs.sniffles-bionano.hmp.txt | grep -v "#" |  grep -v scaf | wc -l`
echo "sv total $svcount"

svcds=`cat data/variants/NAM-structural-variations_CDS.txt | grep -v "#" |  grep -v scaf | wc -l`
echo "sv total $svcds"

awk '$5 > $4' data/ref/zea_maysb73_core_3_87_1_chr-added.gff > TMP
snpcds=`awk '$3 ~ /CDS/' TMP | grep -v scaf | bedtools intersect -wa -a data/variants/B73v5.NAM-illumina_filtered-pass-only-two-round-gatk-snps.vcf.gz -b stdin |  grep -v scaf | wc -l`
echo "number of cds snps $snpcds"

rm TMP

