#!/bin/bash
ref=$1
transcript=$2
nam_ref=${ref#ref_genomes/}
nam_tran=${transcript#transcript_fasta/}
nam=${nam_ref%.}_${nam_tran%.fasta}
module load blast-plus/2.7.1-py2-vvbzyor #loads the module

blastn -subject $ref -query $transcript -out out_${nam} -outfmt "6 qseqid sseqid pident qstart qend length sstart send qcovs" 
#blastn command where ref = subject, transcript fasta = query, name of output file is combo of both
#6 = tabular, #qseqid = query seq ID, #sseqid = subject seq ID
#pident = percentage identical matches 
#qstart start of alignment in query, #qend end of alignment in query (S = subject for sstart and send)
#length = Alignment length 
#qcovs= query coverage per subject