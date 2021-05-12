#!/bin/bash
module load star
cd /work/LAS/mhufford-lab/arnstrm/newNAM/analyses/f-rnaseq/index
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/plant/Zea_mays/latest_assembly_versions/GCA_000005005.6_B73_RefGen_v4/GCA_000005005.6_B73_RefGen_v4_genomic.fna.gz
gunzip GCA_000005005.6_B73_RefGen_v4_genomic.fna.gz
genome=GCA_000005005.6_B73_RefGen_v4_genomic.fna
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/plant/Zea_mays/latest_assembly_versions/GCA_000005005.6_B73_RefGen_v4/GCA_000005005.6_B73_RefGen_v4_genomic.gff.gz
gunzip GCA_000005005.6_B73_RefGen_v4_genomic.gff.gz
GFF=GCA_000005005.6_B73_RefGen_v4_genomic.gff
index="/work/LAS/mhufford-lab/arnstrm/newNAM/analyses/f-rnaseq/index/B73.v4"

mkdir -p $index

STAR \
  --runMode genomeGenerate \
  --runThreadN 36 \
  --genomeDir $index \
  --genomeFastaFiles ${genome}
