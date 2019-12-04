#!/bin/bash
if [ "$#" -ne 2 ] ; then
echo "please provide:"
echo -e "\t\t(1) name for index, preferably the NAM line name"
echo -e "\t\t(2) genome sequences, use only scaffolds"
echo "";
echo "./00_build_index_hisat2 <NAM name> <FASTA-scaffolds-file>" ;
echo "";
exit 0;
fi

NAM=$(basename $(dirname $(pwd)))
file=$2

module load hisat2
hisat2-build ${file} $NAM
