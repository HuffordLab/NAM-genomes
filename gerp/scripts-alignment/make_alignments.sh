ref=./data/sequences/ref/Zm-B73-REFERENCE-NAM-5.0.fa
ref_name=$( basename $ref )
db_prefix=./analyses/last/lastdb/${ref_name}-MAM4
query=$1
query_name=$(basename $query .gz)
query_name=$(basename $query_name .fa)

mkdir -p ./analyses/last/mat/
mat=./analyses/last/mat/"$ref_name"_"$query_name".mat

last-train -P0 --revsym --matsym --gapsym -E0.05 -C2 $db_prefix $query > $mat

mkdir -p ./analyses/last/maf/

maf=./analyses/last/maf/Zm-B73-REFERENCE-NAM-5.0_"$query_name".maf

lastal -m50 -E0.05 -C2 -p $mat $db_prefix $query > $maf
# m50 makes allowed multiplicity of initial hits 50, each match lengthened
# until it occurs at most this many times
# E is threshold for errors
# C2 makes lastal discard any gapless alignment in 2 or more gapped alignments
# p specifies match/mismatch score matrix
# Then reference db
# And finally the query fasta


mkdir -p ./analyses/last/axt
axt=./analyses/last/axt/"$ref_name"_"$query_name".axt

maf-convert axt $maf > $axt

mkdir -p ./analyses/last/chain
chain=./analyses/last/chain/"$ref_name"_"$query_name".chain
merged_chain=./analyses/last/chain/"$ref_name"_"$query_name".all.chain

axtChain $axt $ref $query $chain -linearGap=loose -faQ -faT

# axtChain options
# linearGap=loose for distant species
# -faQ The specified qNibDir is a fasta file with multiple sequences for query
# -faT The specified tNibDir is a fasta file with multiple sequences for target

mkdir -p ./analyses/last/chain_merged
merged_chain=./analyses/last/chain/"$ref_name"_"$query_name".all.chain

chainMergeSort $chain > $merged_chain

mkdir -p ./analyses/chromsize/
ref_size=./analyses/chromsize/"$ref_name".size
query_size=./analyses/chromsize/"$query_name".size

if [ ! -f $ref_size ]; then
  faSize $ref -detailed > $ref_size
fi

if [ ! -f $query_size ]; then
  faSize $query -detailed > $query_size
fi

# note: if your chromosomes are id'd with numbers, faSize will add everything
# before .fa in the file name to the beginning of each chromosome id in the size
# file. This causes issues in the next step. Either make all id's non-numeric
# (e.g. 'chr1' instead of '1') or go in to size file and remove prefix in front
# of chromosome ids.

mkdir -p ./analyses/last/chain_prenet/
chain_prenet=./analyses/last/chain_prenet/"$ref_name"_"$query_name".all.pre.chain
chainPreNet $merged_chain $ref_size $query_size $chain_prenet

mkdir -p ./analyses/last/target_net/
mkdir -p ./analyses/last/query_net/
target_net=./analyses/last/target_net/"$ref_name"_"$query_name".net
query_net=./analyses/last/query_net/"$query_name"_"$ref_name".net
chainNet $chain_prenet $ref_size $query_size $target_net $query_net

# making 2bit files for netToAxt
mkdir -p ./analyses/last/2bit/
query_twobit=./analyses/last/2bit/"$query_name".2bit
ref_twobit=./analyses/last/2bit/"$ref_name".2bit

if [ ! -f $query_twobit ]; then
  faToTwoBit $query $query_twobit
fi

if [ ! -f $ref_twobit ]; then
  faToTwoBit $ref $ref_twobit
fi

mkdir -p ./analyses/last/net_axt/
net_axt=./analyses/last/net_axt/"$ref_name"_"$query_name".net.axt
netToAxt $target_net $chain_prenet $ref_twobit $query_twobit $net_axt

mkdir -p ./analyses/last/net_axt/net_maf
net_maf=./analyses/last/net_axt/net_maf/"$ref_name"_"$query_name".net.maf

axtToMaf $net_axt $ref_size $query_size $net_maf

# now to make mafs one to one

head -n 29 $maf > "$net_maf"_w_header
cat $net_maf >> "$net_maf"_w_header

one_maf=./analyses/last/net_axt/net_maf/"$ref_name"_"$query_name".1to1.maf

last-split -m1 "$net_maf"_w_header |
maf-swap |
awk -v q="$query_name" -v r="$ref_name" '/^s/ {$2 = (++s % 2 ? q "." : r ".") \
$2} 1' | last-split -m1 | \
maf-swap | last-postmask > $one_maf
