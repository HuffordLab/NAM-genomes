#!/bin/bash
module load BEDTools/2.28.0-foss-2018a

gap_dir=/scratch/jl03308/NAM_Canu1.8_verified_version_1/pseudomolecules/Ngap
TE_dir=/scratch/jl03308/NAM_Canu1.8_verified_version_1/annotation/pan_TE_v1.0
gene_dir=/scratch/jl03308/NAM_Canu1.8_verified_version_1/annotation/NAM_genome_and_annotation_Jan2020_release
wdir=/scratch/jl03308/NAM_pancentromere/analysis/peak_call/NAM

cd ${TE_dir}
for file in *.fasta.mod.EDTA.TEanno.gff; do line=$(basename $file |cut -f1 -d ".")
cat $file |awk -v srch="${line}_" -v repl="" '{ sub(srch,repl,$0); print $0 }'|awk '{if($9~"Gypsy_huck"){print$0}}' |cut -f1,4,5 > ${line}.huck.bed
cat $file |awk -v srch="${line}_" -v repl="" '{ sub(srch,repl,$0); print $0 }'|awk '{if($9~"Gypsy_cinful_zeon"){print$0}}' |cut -f1,4,5 > ${line}.cinful-zeon.bed
cat $file |awk -v srch="${line}_" -v repl="" '{ sub(srch,repl,$0); print $0 }'|awk '{if($9~"Gypsy_grande"){print$0}}' |cut -f1,4,5 > ${line}.grande.bed
cat $file |awk -v srch="${line}_" -v repl="" '{ sub(srch,repl,$0); print $0 }'|awk '{if($9~"Gypsy_prem1"||$9~"Gypsy_xilon_diguus"||$9~"Gypsy_tekay"){print$0}}' |cut -f1,4,5 > ${line}.prem1.bed
cat $file |awk -v srch="${line}_" -v repl="" '{ sub(srch,repl,$0); print $0 }'|awk '{if($9~"Copia_opie"||$9~"Copia_ji"||$9~"Copia_ruda"||$9~"Copia_giepum"){print$0}}' |cut -f1,4,5 > ${line}.opie-ji.bed
cat $file |awk -v srch="${line}_" -v repl="" '{ sub(srch,repl,$0); print $0 }' | awk '{if($9~"CRM1"){print$0}}' |cut -f1,4,5 > $line.crm1.bed
cat $file |awk -v srch="${line}_" -v repl="" '{ sub(srch,repl,$0); print $0 }' | awk '{if($9~"CRM2"){print$0}}' |cut -f1,4,5 > $line.crm2.bed
cat $file | awk -v srch="${line}_" -v repl="" '{ sub(srch,repl,$0); print $0 }' | awk '{if($3~"DNA"||$3~"LINE"||$3~"LTR"||$3~"MITE"){print$0}}'|cut -f1,4,5 > $line.allTE.bed
cat $file|awk -v srch="${line}_" -v repl="" '{ sub(srch,repl,$0); print $0 }' | awk '{if($9~"Gypsy_flip"){print$0}}' |cut -f1,4,5 > $line.flip.bed
cat $file|awk -v srch="${line}_" -v repl="" '{ sub(srch,repl,$0); print $0 }' | awk '{if($9~"Gypsy_gyma"){print$0}}' |cut -f1,4,5 > $line.gyma.bed
cat $file|awk -v srch="${line}_" -v repl="" '{ sub(srch,repl,$0); print $0 }' | awk '{if($9~"Gypsy_uwum"){print$0}}' |cut -f1,4,5 > $line.uwum.bed
cat $file|awk -v srch="${line}_" -v repl="" '{ sub(srch,repl,$0); print $0 }' | awk '{if($9~"Gypsy_leviathan"){print$0}}' |cut -f1,4,5 > $line.leviathan.bed
cat $file|awk -v srch="${line}_" -v repl="" '{ sub(srch,repl,$0); print $0 }' | awk '{if($9~"CRM3"){print$0}}' |cut -f1,4,5 > $line.crm3.bed
cat $file|awk -v srch="${line}_" -v repl="" '{ sub(srch,repl,$0); print $0 }' | awk '{if($9~"ansuya"){print$0}}' |cut -f1,4,5 > $line.ansuya.bed
cat $file|awk -v srch="${line}_" -v repl="" '{ sub(srch,repl,$0); print $0 }' | awk '{if($9~"eugene"){print$0}}' |cut -f1,4,5 > $line.eugene.bed
cat $file |awk -v srch="${line}_" -v repl="" '{ sub(srch,repl,$0); print $0 }'|awk '{if($9~"doke"){print$0}}' |cut -f1,4,5 > ${line}.doke.bed
cat $file |awk -v srch="${line}_" -v repl="" '{ sub(srch,repl,$0); print $0 }'|awk '{if($9~"Copia"){print$0}}' |cut -f1,4,5 > ${line}.Copia.bed
cat $file |awk -v srch="${line}_" -v repl="" '{ sub(srch,repl,$0); print $0 }'|awk '{if($9~"iteki"){print$0}}' |cut -f1,4,5 > ${line}.iteki.bed
cat $file|awk -v srch="${line}_" -v repl="" '{ sub(srch,repl,$0); print $0 }' | awk '{if($9~"CRM" || $9~"ansuya"){print$0}}' |cut -f1,4,5 |bedtools sort -i - |bedtools merge -i - > $line.crm_ansuya.bed

#cat ${line}.crm.bed $line.ansuya.bed |bedtools sort -i - |bedtools merge -i - > $line.crm_ansuya.bed
done

#cat *.pseudomolecules-v1.fasta.mod.EDTA.TEanno.chr.gff > NAM.pseudomolecules-v1.fasta.mod.EDTA.TEanno.chr.gff


cd ${gene_dir}
for file in */Zm*.gff; do line=$(basename $file |cut -f2 -d "-")
cat $file | grep -v "#"|awk '{if($3=="gene"){print$0}}' |cut -f1,4,5  > $line.evd.bed
done


cd ${wdir}
lines=$(cat NAM.centro.mean.coords |cut -f1 |sort |uniq)

'''map gap info to centromere range'''
for line in $lines 
do
cat NAM.centro.mean.coords |awk -v var1=$line '{if($1==var1){print$0}}'  |cut -f2,3,4,5 |bedtools intersect -c -a - -b ${gap_dir}/${line}.100n.bed | bedtools intersect -c -a -  -b ${gap_dir}/${line}.13n.bed| bedtools intersect -wao -a -  -b ${gap_dir}/${line}.norm_n.bed|cut -f1,2,3,4,5,6,11 | awk -v var1=$line '{print var1"\t"$0}'|sort -k1,1  >> NAM.centro.coords.gaps
done

'''check assembly stats by size (not assembled: smaller than 1500000 or contains a 100N gap)'''
cat NAM.centro.coords.gaps |awk '{if($5<1500000||$6>0){print$0"\tN"} else{print$0"\tY"}}' |bedtools groupby -i - -g 1,2,3,4,5,6,7,9 -c 8 -o sum |awk -v OFS='\t' '{print$1,$2,$3,$4,$5,$6,$7,$9,$8}' > NAM.centro.coords.assembly.status

'''obtain all location of known sequences (centc, gaps, TE and genes)'''
for line in $lines
do
cat ${TE_dir}/${line}.centC.bed ${TE_dir}/${line}.allTE.bed ${gene_dir}/${line}.evd.bed ${gap_dir}/${line}.ngap.bed | cut -f1,2,3 |sort -k1,1 -k2,2n | bedtools merge -i - > ${line}.known_sequences.bed
done

'''check centC, crm2, crm1, and the five major TE content in centromere'''
# for line in $lines
# do
# cat NAM.centro.coords.assembly.status |awk -v var1=$line '{if($1==var1){print$0}}'  |cut -f2,3,4,5,6,7,8,9 | \
# bedtools intersect -a - -b ${TE_dir}/${line}.centC.bed -wao| bedtools groupby -i - -g 1,2,3,4,5,6,7,8 -c 12 -o sum | \
# bedtools intersect -a - -b ${TE_dir}/${line}.crm2.bed -wao| bedtools groupby -i - -g 1,2,3,4,5,6,7,8,9 -c 13 -o sum | \
# bedtools intersect -a - -b ${TE_dir}/${line}.crm1.bed -wao| bedtools groupby -i - -g 1,2,3,4,5,6,7,8,9,10 -c 14 -o sum | \
# bedtools intersect -a - -b ${TE_dir}/${line}.cinful-zeon.bed -wao| bedtools groupby -i - -g 1,2,3,4,5,6,7,8,9,10,11 -c 15 -o sum | \
# bedtools intersect -a - -b ${TE_dir}/${line}.opie-ji.bed -wao| bedtools groupby -i - -g 1,2,3,4,5,6,7,8,9,10,11,12 -c 16 -o sum | \
# bedtools intersect -a - -b ${TE_dir}/${line}.huck.bed -wao| bedtools groupby -i - -g 1,2,3,4,5,6,7,8,9,10,11,12,13 -c 17 -o sum | \
# bedtools intersect -a - -b ${TE_dir}/${line}.prem1.bed -wao| bedtools groupby -i - -g 1,2,3,4,5,6,7,8,9,10,11,12,13,14 -c 18 -o sum | \
# bedtools intersect -a - -b ${TE_dir}/${line}.grande.bed -wao| bedtools groupby -i - -g 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15 -c 19 -o sum | \
# bedtools intersect -a - -b ${TE_dir}/${line}.allTE.bed -wao| bedtools groupby -i - -g 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16 -c 20 -o sum | \
# bedtools intersect -a - -b ${gene_dir}/${line}.evd.bed -wao| bedtools groupby -i - -g 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17 -c 21 -o sum | \
# bedtools intersect -a - -b ${line}.known_sequences.bed -wao| bedtools groupby -i - -g 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18 -c 22 -o sum | \
# awk -v var1=$line '{print var1"\t"$0}'|sort -k1,1 -k2,2 >> NAM.centro.coords.assembly.status.content
# done
# #1) other TE  (allTEs-crm2 -crm1 -5 major families) 
# #2) total gap (100, 13, known)
# #3) unknown (total_len-known)
# 
# cat NAM.centro.coords.assembly.status.content| awk '{print$0"\t"$18-$11-$12-$13-$14-$15-$16-$17"\t"$6*100+$8+$7*13"\t"$5-$20}' > NAM.centro.coords.assembly.status.content.sum
# 
for line in $lines
do
cat NAM.centro.coords.assembly.status |awk -v var1=$line '{if($1==var1){print$0}}'  |cut -f2,3,4,5,6,7,8,9 | \
bedtools intersect -a - -b ${TE_dir}/${line}.centC.bed -wao| bedtools groupby -i - -g 1,2,3,4,5,6,7,8 -c 12 -o sum | \
bedtools intersect -a - -b ${TE_dir}/${line}.crm_ansuya.bed -wao| bedtools groupby -i - -g 1,2,3,4,5,6,7,8,9 -c 13 -o sum | \
bedtools intersect -a - -b ${TE_dir}/${line}.flip.bed -wao| bedtools groupby -i - -g 1,2,3,4,5,6,7,8,9,10 -c 14 -o sum | \
bedtools intersect -a - -b ${TE_dir}/${line}.prem1.bed -wao| bedtools groupby -i - -g 1,2,3,4,5,6,7,8,9,10,11 -c 15 -o sum | \
bedtools intersect -a - -b ${TE_dir}/${line}.cinful-zeon.bed -wao| bedtools groupby -i - -g 1,2,3,4,5,6,7,8,9,10,11,12 -c 16 -o sum | \
bedtools intersect -a - -b ${TE_dir}/${line}.gyma.bed -wao| bedtools groupby -i - -g 1,2,3,4,5,6,7,8,9,10,11,12,13 -c 17 -o sum | \
bedtools intersect -a - -b ${TE_dir}/${line}.doke.bed -wao| bedtools groupby -i - -g 1,2,3,4,5,6,7,8,9,10,11,12,13,14 -c 18 -o sum | \
bedtools intersect -a - -b ${TE_dir}/${line}.iteki.bed -wao| bedtools groupby -i - -g 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15 -c 19 -o sum | \
bedtools intersect -a - -b ${TE_dir}/${line}.grande.bed -wao| bedtools groupby -i - -g 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16 -c 20 -o sum | \
bedtools intersect -a - -b ${TE_dir}/${line}.leviathan.bed -wao| bedtools groupby -i - -g 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17 -c 21 -o sum | \
bedtools intersect -a - -b ${TE_dir}/${line}.eugene.bed -wao| bedtools groupby -i - -g 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18 -c 22 -o sum | \
bedtools intersect -a - -b ${TE_dir}/${line}.uwum.bed -wao| bedtools groupby -i - -g 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19 -c 23 -o sum | \
bedtools intersect -a - -b ${TE_dir}/${line}.Copia.bed -wao| bedtools groupby -i - -g 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20 -c 24 -o sum | \
bedtools intersect -a - -b ${TE_dir}/${line}.allTE.bed -wao| bedtools groupby -i - -g 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21 -c 25 -o sum | \
bedtools intersect -a - -b ${gene_dir}/${line}.evd.bed -wao| bedtools groupby -i - -g 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22 -c 26 -o sum | \
bedtools intersect -a - -b ${line}.known_sequences.bed -wao| bedtools groupby -i - -g 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23 -c 27 -o sum | \
awk -v var1=$line '{print var1"\t"$0}'|sort -k1,1 -k2,2 >> NAM.centro.coords.assembly.status.content
done

#$10-$25:centc, crm, flip,prem1, cinful-zeon, gyma, doke, iteki,grande, leviathan, centA, uwum, Copia, all TE, gene, known_sequences
#1) other TE  (allTEs-crm ... copia) 
#2) total gap (100, 13, known)
#3) unknown (total_len-known)

cat NAM.centro.coords.assembly.status.content| awk '{print$0"\t"$23-$11-$12-$13-$14-$15-16-$17-$18-$19-$20-$21-$22"\t"$6*100+$8+$7*13"\t"$5-$25}' > NAM.centro.coords.assembly.status.content.sum
