#!/bin/bash
#This script will calculate the assembly stats and BUSCO metrics for the input genome

python3 assembly_stats.py $1 $2 
python run_BUSCO.py -i $1 -l $4 -o $3 -m genome -sp $5 -c $6 --blast_single_core 
