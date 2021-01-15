#!/bin/bash
window="/work/LAS/mhufford-lab/arnstrm/newNAM/analyses/l-gatk-illumina/b-genome/B73.PLATINUM.pseudomolecules-v1_1mb_coords.bed"

while read line; do
out=$(echo $line |sed 's/:/_/g')
echo "java -Djava.io.tmpdir=\$TMPDIR -Xmx20G -jar /work/LAS/mhufford-lab/arnstrm/programs/gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar HaplotypeCaller -R /work/LAS/mhufford-lab/arnstrm/newNAM/analyses/l-gatk-illumina/b-genome/B73.PLATINUM.pseudomolecules-v1.fasta -I 1487-1_S4_L001_R1_001.fastq_final_recal.bam -I 1487-1_S4_L002_R1_001.fastq_final_recal.bam -I 1487-1_S4_L003_R1_001.fastq_final_recal.bam -I 1487-1_S4_L004_R1_001.fastq_final_recal.bam -I 1487-2_S2_L001_R1_001.fastq_final_recal.bam -I 1487-2_S2_L002_R1_001.fastq_final_recal.bam -I 1487-2_S2_L003_R1_001.fastq_final_recal.bam -I 1487-2_S2_L004_R1_001.fastq_final_recal.bam -I 1487-3_S1_L001_R1_001.fastq_final_recal.bam -I 1487-3_S1_L002_R1_001.fastq_final_recal.bam -I 1487-3_S1_L003_R1_001.fastq_final_recal.bam -I 1487-3_S1_L004_R1_001.fastq_final_recal.bam -I 1487-4_S3_L001_R1_001.fastq_final_recal.bam -I 1487-4_S3_L002_R1_001.fastq_final_recal.bam -I 1487-4_S3_L003_R1_001.fastq_final_recal.bam -I 1487-4_S3_L004_R1_001.fastq_final_recal.bam -I 1506-1_S1_L001_R1_001.fastq_final_recal.bam -I 1506-1_S1_L002_R1_001.fastq_final_recal.bam -I 1506-1_S1_L003_R1_001.fastq_final_recal.bam -I 1506-1_S1_L004_R1_001.fastq_final_recal.bam -I 1507-1_S1_L001_R1_001.fastq_final_recal.bam -I 1507-1_S1_L002_R1_001.fastq_final_recal.bam -I 1507-1_S1_L003_R1_001.fastq_final_recal.bam -I 1507-1_S1_L004_R1_001.fastq_final_recal.bam -I 1508-1_S1_L001_R1_001.fastq_final_recal.bam -I 1508-1_S1_L002_R1_001.fastq_final_recal.bam -I 1508-1_S1_L003_R1_001.fastq_final_recal.bam -I 1508-1_S1_L004_R1_001.fastq_final_recal.bam -I 1605-1_S1_L001_R1_001.fastq_final_recal.bam -I 1605-1_S1_L002_R1_001.fastq_final_recal.bam -I 1605-1_S1_L003_R1_001.fastq_final_recal.bam -I 1605-1_S1_L004_R1_001.fastq_final_recal.bam -I 1721-1_S1_L001_R1_001.fastq_final_recal.bam -I 1721-1_S1_L002_R1_001.fastq_final_recal.bam -I 1721-1_S1_L003_R1_001.fastq_final_recal.bam -I 1721-1_S1_L004_R1_001.fastq_final_recal.bam -I 1721-2_S1_L001_R1_001.fastq_final_recal.bam -I 1721-2_S1_L002_R1_001.fastq_final_recal.bam -I 1721-2_S1_L003_R1_001.fastq_final_recal.bam -I 1721-2_S1_L004_R1_001.fastq_final_recal.bam -I 1721-5_S1_L001_R1_001.fastq_final_recal.bam -I 1721-5_S1_L002_R1_001.fastq_final_recal.bam -I 1721-5_S1_L003_R1_001.fastq_final_recal.bam -I 1721-5_S1_L004_R1_001.fastq_final_recal.bam -I 1726-Oh7B_S1_L001_R1_001.fastq_final_recal.bam -I 1726-Oh7B_S1_L002_R1_001.fastq_final_recal.bam -I 1726-Oh7B_S1_L003_R1_001.fastq_final_recal.bam -I 1726-Oh7B_S1_L004_R1_001.fastq_final_recal.bam -I 2025-B73_S1_L001_R1_001.fastq_final_recal.bam -I 2025-B73_S1_L002_R1_001.fastq_final_recal.bam -I 2025-B73_S1_L003_R1_001.fastq_final_recal.bam -I 2025-B73_S1_L004_R1_001.fastq_final_recal.bam -I 2025-NC350_S1_L001_R1_001.fastq_final_recal.bam -I 2025-NC350_S1_L002_R1_001.fastq_final_recal.bam -I 2025-NC350_S1_L003_R1_001.fastq_final_recal.bam -I 2025-NC350_S1_L004_R1_001.fastq_final_recal.bam -I B73Ab10_R1_part_aa_final_recal.bam -I B73Ab10_R1_part_ab_final_recal.bam -I B73Ab10_R1_part_ac_final_recal.bam -I B73Ab10_R1_part_ad_recal.bam -I B73Ab10_R1_part_ae_recal.bam -I B73Ab10_R1_part_af_recal.bam -I B73Ab10_R1_part_ag_recal.bam -I CML103-1_S1_L001_R1_001.fastq_final_recal.bam -I CML103-1_S1_L002_R1_001.fastq_final_recal.bam -I CML103-1_S1_L003_R1_001.fastq_final_recal.bam -I CML103-1_S1_L004_R1_001.fastq_final_recal.bam -I CML333_S0_L001_R1_001.fastq_final_recal.bam -I CML333_S0_L002_R1_001.fastq_final_recal.bam -I CML333_S0_L003_R1_001.fastq_final_recal.bam -I CML333_S0_L004_R1_001.fastq_final_recal.bam -I HP301_R1_part_aa_final_recal.bam -I HP301_R1_part_ab_final_recal.bam -I HP301_R1_part_ac_final_recal.bam -I HP301_R1_part_ad_final_recal.bam -I HP301_R1_part_ae_final_recal.bam -I Ki11_R1_part_aa_final_recal.bam -I Ki11_R1_part_ab_final_recal.bam -I Ki11_R1_part_ac_final_recal.bam -I Ki11_R1_part_ad_final_recal.bam -I Ki11_R1_part_ae_final_recal.bam -I Ki3-1626_S1_L001_R1_001.fastq_final_recal.bam -I Ki3-1626_S1_L002_R1_001.fastq_final_recal.bam -I Ki3-1626_S1_L003_R1_001.fastq_final_recal.bam -I Ki3-1626_S1_L004_R1_001.fastq_final_recal.bam -I M162W_R1_part_aa_final_recal.bam -I M162W_R1_part_ab_recal.bam -I M162W_R1_part_ac_recal.bam -I M162W_R1_part_ad_recal.bam -I M162W_R1_part_ae_recal.bam -I Maize-10X-1_S4_L001_R1_001.fastq_final_recal.bam -I Maize-10X-1_S4_L002_R1_001.fastq_final_recal.bam -I Maize-10X-1_S4_L003_R1_001.fastq_final_recal.bam -I Maize-10X-1_S4_L004_R1_001.fastq_final_recal.bam -I Maize-10X-2_S2_L001_R1_001.fastq_final_recal.bam -I Maize-10X-2_S2_L002_R1_001.fastq_final_recal.bam -I Maize-10X-2_S2_L003_R1_001.fastq_final_recal.bam -I Maize-10X-2_S2_L004_R1_001.fastq_final_recal.bam -I Maize-10X-3_S3_L001_R1_001.fastq_final_recal.bam -I Maize-10X-3_S3_L002_R1_001.fastq_final_recal.bam -I Maize-10X-3_S3_L003_R1_001.fastq_final_recal.bam -I Maize-10X-3_S3_L004_R1_001.fastq_final_recal.bam -I Maize-10X-4_S1_L001_R1_001.fastq_final_recal.bam -I Maize-10X-4_S1_L002_R1_001.fastq_final_recal.bam -I Maize-10X-4_S1_L003_R1_001.fastq_final_recal.bam -I Maize-10X-4_S1_L004_R1_001.fastq_final_recal.bam -I Mo18W_R1_part_aa_final_recal.bam -I Mo18W_R1_part_ab_final_recal.bam -I Mo18W_R1_part_ac_final_recal.bam -I Mo18W_R1_part_ad_final_recal.bam -I MS37W_R1_part_aa_recal.bam -I MS37W_R1_part_ab_recal.bam -I MS37W_R1_part_ac_recal.bam -I MS37W_R1_part_ad_recal.bam -I MS37W_R1_part_ae_recal.bam -I MS71-1797_S1_L001_R1_001.fastq_final_recal.bam -I MS71-1797_S1_L002_R1_001.fastq_final_recal.bam -I MS71-1797_S1_L003_R1_001.fastq_final_recal.bam -I MS71-1797_S1_L004_R1_001.fastq_final_recal.bam -I NC358_R1_part_aa_final_recal.bam -I NC358_R1_part_ab_final_recal.bam -I NC358_R1_part_ac_final_recal.bam -I OH43-1721_S1_L001_R1_001.fastq_final_recal.bam -I OH43-1721_S1_L002_R1_001.fastq_final_recal.bam -I OH43-1721_S1_L003_R1_001.fastq_final_recal.bam -I OH43-1721_S1_L004_R1_001.fastq_final_recal.bam -I P39-1927_S1_L001_R1_001.fastq_final_recal.bam -I P39-1927_S1_L002_R1_001.fastq_final_recal.bam -I P39-1927_S1_L003_R1_001.fastq_final_recal.bam -I P39-1927_S1_L004_R1_001.fastq_final_recal.bam -I Tx303-1927_S1_L001_R1_001.fastq_final_recal.bam -I Tx303-1927_S1_L002_R1_001.fastq_final_recal.bam -I Tx303-1927_S1_L003_R1_001.fastq_final_recal.bam -I Tx303-1927_S1_L004_R1_001.fastq_final_recal.bam -I Tzi8_S1_L001_R1_001.fastq_final_recal.bam -I Tzi8_S1_L002_R1_001.fastq_final_recal.bam -I Tzi8_S1_L003_R1_001.fastq_final_recal.bam -I Tzi8_S1_L004_R1_001.fastq_final_recal.bam -I Undetermined_S0_L001_R1_001.fastq_final_recal.bam -I Undetermined_S0_L002_R1_001.fastq_final_recal.bam -I Undetermined_S0_L003_R1_001.fastq_final_recal.bam -I Undetermined_S0_L004_R1_001.fastq_final_recal.bam -L $line  -O $out.vcf"
done<${window}