ref_rm=./data/sequences/ref/Zm-B73-REFERENCE-NAM-5.0_rm.fa

ref_rm_name=$( basename $ref_rm )

combined_maf=./analyses/last/net_axt/net_maf/combined.maf

# splitting combined maf by target sequence
mkdir -p ./analyses/last/split_maf/
outroot=./analyses/last/split_maf/
mafSplit -byTarget dummy.bed $outroot "$combined_maf".filtered -useFullSequenceName

# mafSplit throws an error if you don't put dummy.bed, even though byTarget
# means it's being ignored

path_to_phast=~/bin/phast
msa_view=${path_to_phast}/bin/msa_view
maf_dir=($( ls -d ./analyses/last/split_maf/*chr*.maf ))


mkdir -p ./analyses/last/msa_fasta/
mkdir -p ./data/sequences/ref/split/

faSplit byname $ref_rm ./data/sequences/ref/split/

path_to_match_masking=~/projects/nam/scripts/matchMasking.pl

for maf_file in "${maf_dir[@]}"
do
  chr=$(basename $maf_file | sed -e 's/0\(.*\).maf/\1/')

  fasta=./analyses/last/msa_fasta/"$chr".fa
  ref_rm_chr=./data/sequences/ref/split/"$chr".fa
  rm_fasta=./analyses/last/msa_fasta/"$chr"_rm.fa
  $msa_view $maf_file -f -G 1 --refseq $ref_rm_chr > $fasta
  sed -E 's/> />/g' $fasta > "$fasta"_tmp && mv "$fasta"_tmp $fasta
  perl $path_to_match_masking \
  --ref $ref_rm_chr \
  --fasta $fasta \
  --out $rm_fasta
done
