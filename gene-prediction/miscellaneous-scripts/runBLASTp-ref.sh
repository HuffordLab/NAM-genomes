#!/bin/bash
# perfomrs NR blast (blastx)
infile="$1"
outfile="$(basename "${infile%.*}")-refseq.out"
database="/work/GIF/databases/ncbi_refseq_protein/refseq_protein"
module load  blast-plus
blastp \
 -query "${infile}" \
 -db "${database}" \
 -out "${outfile}" \
 -evalue 1e-3 \
 -num_threads 4 \
 -outfmt 6
# -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids salltitles qcovs"
