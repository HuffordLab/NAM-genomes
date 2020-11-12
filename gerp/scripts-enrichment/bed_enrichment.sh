test_bed=./analyses/gerp/all.elems.bed
cat ./analyses/gerp/*.elems.bed >> $test_bed

SV=./data/variants/NAM_founders_SVs.bed
SNP=./data/variants/B73v5.NAM-illumina-snps-only_filtered.frq.bed
genome=./data/chromsize/Zm-B73-REFERENCE-NAM-5.0.size.bed

# working with v5 gff
# get genes
gff=./data/annotations/zea_maysb73_core_3_87_1.gff
sed -i -e 's/^/chr/' $gff
genes=./data/annotations/genes_B73v5.gff3
genes_merged=./data/annotations/genes_B73v5_merged.gff3

awk 'BEGIN {OFS="\t"} $3~/gene/ {print $0}' $gff |
awk 'BEGIN {OFS="\t"} $5>$4 {print $0}'> $genes
sort -k 1,1 -k4,4n $genes > "$genes"_tmp && mv "$genes"_tmp $genes
bedtools merge -i $genes > $genes_merged
cds=/home/aihudson/projects/nam/data/annotations/cds_B73v5.gff3
awk 'BEGIN {OFS="\t"} $3=="CDS" {print $0}' $gff |
awk 'BEGIN {OFS="\t"} $5>$4 {print $0}' > $cds

# coding - get bed file with just coding regions

gerp_cds=./analyses/gerp/gerp_cds.elems.bed
bedtools intersect -a $test_bed -b $cds > $gerp_cds
sort -k 1,1 -k2,2n $gerp_cds > "$gerp_cds"_tmp && mv "$gerp_cds"_tmp $gerp_cds
bedtools merge -i $gerp_cds > "$gerp_cds"_tmp && mv "$gerp_cds"_tmp $gerp_cds

# get non-genic regions
gerp_nongenic=./analyses/gerp/gerp_nongenic.elems.bed
bedtools intersect -v -a $test_bed -b $genes_merged > $gerp_nongenic
sort -k 1,1 -k2,2n $gerp_nongenic > "$gerp_nongenic"_tmp && mv "$gerp_nongenic"_tmp $gerp_nongenic
bedtools merge -i $gerp_nongenic > "$gerp_nongenic"_tmp && mv "$gerp_nongenic"_tmp $gerp_nongenic

# get non-coding regions
gerp_noncds=./analyses/gerp/gerp_noncds.elems.bed
bedtools intersect -v -a $test_bed -b $cds > $gerp_noncds
sort -k 1,1 -k2,2n $gerp_noncds > "$gerp_noncds"_tmp && mv "$gerp_noncds"_tmp $gerp_noncds
bedtools merge -i $gerp_noncds > "$gerp_noncds"_tmp && mv "$gerp_noncds"_tmp $gerp_noncds

sort -k 1,1 -k2,2n $SV > "$SV"_tmp && mv "$SV"_tmp $SV

SV_merged=./data/variants/NAM_founders_SVs_merged.bed
sort -k 1,1 -k2,2n $SNP > "$SNP"_tmp && mv "$SNP"_tmp $SNP
SNP_merged=./data/variants/B73v5.NAM-illumina-snps-only_filtered_merged.frq.bed
bedtools merge -i $SNP > $SNP_merged

sh_helper_subset=./scripts/scripts-enrichment/bed_enrichment_helper_subset.sh
sh_helper=./scripts/scripts-enrichment/bed_enrichment_helper.sh
r_helper=./scripts/scripts-enrichment/bed_enrichment_helper.R

# Are structural variants enriched in the bed file?
chi=./analyses/chi/gerp_sv_chi.txt
bash $sh_helper $test_bed $SV_merged $genome $chi
# #
# #
# #
Rscript $r_helper $chi

# Are deletions enriched in the bed file?
del=./data/variants/NAM_founders_del.bed
sort -k 1,1 -k2,2n $del > "$del"_tmp && mv "$del"_tmp $del
del_merged=./data/variants/NAM_founders_del_merged.bed
bedtools merge -i $del > $del_merged
chi=./analyses/chi/gerp_del_chi.txt
bash $sh_helper $test_bed $del_merged $genome $chi
Rscript $r_helper $chi

# Are insertions enriched in the bed file?
ins=./data/variants/NAM_founders_ins.bed
sort -k 1,1 -k2,2n $ins > "$ins"_tmp && mv "$ins"_tmp $ins
ins_merged=./data/variants/NAM_founders_ins_merged.bed
bedtools merge -i $ins > $ins_merged
chi=./analyses/chi/gerp_ne_ins_bp_chi.txt
bash $sh_helper $test_bed $ins_merged $genome $chi
Rscript $r_helper $chi


# Are translocation enriched in the bed file?
tra=./data/variants/NAM_founders_tra.bed
sort -k 1,1 -k2,2n $tra > "$tra"_tmp && mv "$tra"_tmp $tra
tra_merged=./data/variants/NAM_founders_tra_merged.bed
bedtools merge -i $tra > $tra_merged
chi=./analyses/chi/gerp_ne_tra_chi.txt
bash $sh_helper $test_bed $tra_merged $genome $chi
Rscript $r_helper $chi


# Are inversions enriched in the bed file?
inv=./data/variants/NAM_founders_inv.bed
sort -k 1,1 -k2,2n $inv > "$inv"_tmp && mv "$inv"_tmp $inv
inv_merged=./data/variants/NAM_founders_inv_merged.bed
bedtools merge -i $inv > $inv_merged
chi=./analyses/chi/gerp_inv_chi.txt
bash $sh_helper $test_bed $inv_merged $genome $chi
Rscript $r_helper $chi

# Are SNPs enriched in the bed file?

chi=./analyses/chi/gerp_snp_chi.txt
bash $sh_helper $test_bed $SNP $genome $chi
#
#
#
Rscript $r_helper $chi
#
# Are structural variants enriched in coding regions?

chi=./analyses/chi/gerp_cds_sv_chi.txt
gerp_cds=./analyses/gerp/gerp_cds.elems.bed
# bedtools intersect -a $test_bed -b $cds > $gerp_cds
# sort -k 1,1 -k2,2n $gerp_cds > "$gerp_cds"_tmp && mv "$gerp_cds"_tmp $gerp_cds
# bedtools merge -i $gerp_cds > "$gerp_cds"_tmp && mv "$gerp_cds"_tmp $gerp_cds
bash $sh_helper_subset $gerp_cds $SV_merged $genome $chi $test_bed
# #
# #
# #
Rscript $r_helper $chi


# # Are SNPs enriched in coding regions?

chi=./analyses/chi/gerp_cds_snp_chi.txt

bash $sh_helper_subset $gerp_cds $SNP $genome $chi $test_bed

Rscript $r_helper $chi

# Are SVs enriched in non-genic GERP regions?

chi=./analyses/chi/gerp_nongenic_sv_chi.txt

bash $sh_helper_subset $gerp_nongenic $SV_merged $genome $chi $test_bed

Rscript $r_helper $chi

# Are SNPs enriched in non-genic GERP regions?

chi=./analyses/chi/gerp_nongenic_snp_chi.txt
bash $sh_helper_subset $gerp_nongenic $SNP $genome $chi $test_bed

Rscript $r_helper $chi


# Are SVs enriched in non-coding regions?

chi=./analyses/chi/gerp_noncds_sv_chi.txt

bash $sh_helper_subset $gerp_noncds $SV $genome $chi $test_bed



Rscript $r_helper $chi

# Are SNPs enriched in non-coding regions?

chi=./analyses/chi/gerp_noncds_snp_chi.txt

bash $sh_helper_subset $gerp_noncds $SNP $genome $chi $test_bed
#
Rscript $r_helper $chi

mkdir -p ./analyses/windows
windows=./analyses/windows/windows
bedtools makewindows -b $genome -w 10000 > "$windows".bed

bedtools intersect -wao -a $windows -b $ins_merged > "$windows"_ins.bed
bedtools intersect -wao -a $windows -b $del_merged > "$windows"_del.bed
bedtools intersect -wao -a $windows -b $tra_merged > "$windows"_tra.bed
bedtools intersect -wao -a $windows -b $inv_merged > "$windows"_inv.bed
bedtools intersect -wao -a $windows -b $test_bed > "$windows"_gerp.bed
open=./data/annotations/open_chromatin.bed
bedtools intersect -wao -a $windows -b $open > "$windows"_open.bed
masked=./data/annotations/Zm-B73-REFERENCE-NAM-5.0-repeatmask.gff3.gz
masked_merged=./data/annotations/Zm-B73-REFERENCE-NAM-5.0-repeatmask_merged.bed
bedtools merge -i $masked > $masked_merged
bedtools intersect -wao -a $windows -b $masked > "$windows"_masked.bed
