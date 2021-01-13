#!/usr/local/apps/eb/Python/3.7.0-foss-2018a/bin/python

# Import necessary packages
import argparse
import re
from Bio import SeqIO

# Parse command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument("fasta")
args = parser.parse_args()

# Open FASTA, search for masked regions, print in BED3 format
with open(args.fasta) as handle:
    for record in SeqIO.parse(handle, "fasta"):
        for match in re.finditer('N+', str(record.seq)):
            print(record.id, match.start(), match.end(), sep='\t')
