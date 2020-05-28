
# Check if fasta files are zipped, and if they are, unzip them
if [ -f "input_files/hg38_chr1.fa.gz" ] && [ ! -f "input_files/hg38_chr1.fa" ]; then
    echo "Unzipping chr1 fasta file"
    gunzip input_files/hg38_chr1.fa.gz
fi

if [ -f "input_files/hg38_chr11.fa.gz" ] && [ ! -f "input_files/hg38_chr11.fa" ]; then
    echo "Unzipping chr11 fasta file"
    gunzip input_files/hg38_chr11.fa.gz
fi
