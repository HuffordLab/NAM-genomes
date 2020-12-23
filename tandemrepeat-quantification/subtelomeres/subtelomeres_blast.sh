#!/bin/bash

#This script blast the 7 maize subtelomeric sequences to each chromosome of the genome 

module load blast-plus


for file in chr*fasta; do
outname=$(echo $file| cut -d'.' -f 1)
blastn -subject $file -query seq1.fasta -out seq1output_$outname -outfmt "6 qseqid sseqid pident qstart qend length sstart send qcovs"
blastn -subject $file -query seq2.fasta -out seq2output_$outname -outfmt "6 qseqid sseqid pident qstart qend length sstart send qcovs"
blastn -subject $file -query seq3.fasta -out seq3output_$outname -outfmt "6 qseqid sseqid pident qstart qend length sstart send qcovs"
blastn -subject $file -query seq4.fasta -out seq4output_$outname -outfmt "6 qseqid sseqid pident qstart qend length sstart send qcovs"
blastn -subject $file -query seq5.fasta -out seq5output_$outname -outfmt "6 qseqid sseqid pident qstart qend length sstart send qcovs"
blastn -subject $file -query seq6.fasta -out seq6output_$outname -outfmt "6 qseqid sseqid pident qstart qend length sstart send qcovs"
blastn -subject $file -query seq7.fasta -out seq7output_$outname -outfmt "6 qseqid sseqid pident qstart qend length sstart send qcovs"
awk '{ if ($6 >= 330 && $3 >= 80) print $0}' seq1output_$outname > $outname.seq1_filt
awk '{ if ($6 >= 376 && $3 >= 80) print $0}' seq2output_$outname > $outname.seq2_filt
awk '{ if ($6 >= 249 && $3 >=80) print $0}'  seq3output_$outname > $outname.seq3_filt
awk '{ if ($6 >= 305  && $3 >=80) print $0}' seq4output_$outname > $outname.seq4_filt
awk '{ if ($6 >= 627 && $3 >= 80) print $0}' seq5output_$outname > $outname.seq5_filt
awk '{ if ($6 >= 979 && $3 >= 80) print $0}' seq6output_$outname > $outname.seq6_filt
awk '{ if ($6 >= 328 && $3 >= 80) print $0}' seq7output_$outname > $outname.seq7_filt

cat $outname.seq1_filt $outname.seq3_filt $outname.seq3_filt $outname.seq4_filt $outname.seq5_filt $outname.seq6_filt $outname.seq7_filt > $outname.all_filtered_output
awk '{print $2"\t"$7"\t"$8}' $outname.all_filtered_output > $outname.bed
awk '{if($3>$2)print}' $outname.bed | awk '{print $1"\t"$2"\t"$3}' > $outname.f.bed
awk '{if($2>$3)print}' $outname.bed | awk '{print $1"\t"$3"\t"$2}' > $outname.r.bed
cat $outname.f.bed $outname.r.bed > $outname.fr.bed
sort -k2 -n $outname.fr.bed > sorted_$outname.bed
done
