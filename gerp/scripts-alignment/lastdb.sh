ref=./data/sequences/ref/Zm-B73-REFERENCE-NAM-5.0.fa
ref_name=$( basename $ref )
db_prefix=./analyses/last/lastdb/${ref_name}-MAM4

mkdir -p ./analyses/last/lastdb/
lastdb -P0 -uMAM4 -R01 $db_prefix $ref
