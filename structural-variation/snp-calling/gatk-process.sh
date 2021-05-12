#!/bin/bash
## this script merges and processes the files created by split and gather method of GATK
## specifically, this is designed for the first round of SNP calling
## start it in your directory that has split vcf files in in
## 2019/09/19
## Arun Seetharam
## <arnstrm@iastate.edu>
base=second-round
# classify
mkdir -p vcf-files_${base} idx-files_${base}
for f in *.vcf; do
mv $f vcf-files_${base}/;
mv $f.idx idx-files_${base}/;
done
# merge
cd vcf-files_${base}/
for f in *.vcf; do
grep -v "^#" $f > ${f%.*}.nohead;
done
rm body.txt &> /dev/null
cat *.nohead >> body.txt
onevcf=$(ls *.vcf |head -n 1)
grep "^#" ${onevcf} > head.txt
rm ../merged_${base}.vcf &> /dev/null
cat head.txt body.txt >> ../merged_${base}.vcf
rm head.txt body.txt *.nohead
cd ..
# separate SNPs and INDELS and sort them
module load vcftools
vcftools --vcf merged_${base}.vcf --keep-only-indels --recode --recode-INFO-all --out ${base}_indels-only
vcftools --vcf merged_${base}.vcf --remove-indels --recode --recode-INFO-all --out ${base}_snps-only
cat ${base}_snps-only.recode.vcf | vcf-sort -t $TMPDIR -p 36 -c > ${base}_sorted-snps-only.vcf
cat ${base}_indels-only.recode.vcf | vcf-sort -t $TMPDIR -p 36 -c >  ${base}_sorted-indels-only.vcf
# get max-depth value for filtering
grep -v "^#" ${base}_indels-only.recode.vcf | cut -f 8 | grep -oe ";DP=.*" |cut -f 2 -d ";" |cut -f 2 -d "=" > dp_indel.txt
grep -v "^#" ${base}_snps-only.recode.vcf | cut -f 8 | grep -oe ";DP=.*" |cut -f 2 -d ";" |cut -f 2 -d "=" > dp_snps.txt
echo "dp stats for VCF"
cat dp_snps.txt | st
MAXDEPTH=$(cat dp_snps.txt | st | grep -v "^N" |awk '{print ($(NF-1)+$NF)*5}')
echo "dp stats for INDEL"
cat dp_indel.txt | st
ref=/work/LAS/mhufford-lab/arnstrm/newNAM/analyses/l-gatk-illumina/b-genome/B73.PLATINUM.pseudomolecules-v1.fasta
# perform filtering on SNPs
java -Djava.io.tmpdir=$TMPDIR -Xmx20G -jar /work/LAS/mhufford-lab/arnstrm/programs/gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar VariantFiltration \
   --reference $ref \
   --variant ${base}_sorted-snps-only.vcf \
   --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 45.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || DP > ${MAXDEPTH}" \
   --filter-name "FAIL" \
   --output ${base}_sorted-snps-only_filtered.vcf
# perform filtering on INDELS
java -Djava.io.tmpdir=$TMPDIR -Xmx20G -jar /work/LAS/mhufford-lab/arnstrm/programs/gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar VariantFiltration \
   --reference $ref \
   --variant ${base}_sorted-indels-only.vcf \
   --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
   --filter-name "FAIL" \
   --output ${base}_sorted-indels-only_filtered.vcf
grep -v "^#" ${base}_sorted-snps-only_filtered.vcf | awk '$7=="PASS"' > t1
grep "^#" ${base}_sorted-snps-only_filtered.vcf > t0
rm ${base}_sorted-snps_filtered-pass-only.vcf > /dev/null
cat t0 t1 >> ${base}_sorted-snps_filtered-pass-only.vcf

