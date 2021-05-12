#!/bin/bash
input="$1"
nam=$2
#make the actual database
makeblastdb -in $input -dbtype nucl -parse_seqids -title ${nam} # -out ${nam}