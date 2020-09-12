import argparse
import numpy as np
import os
import re
import io
import subprocess

#from parse_read import parse_read

parser = argparse.ArgumentParser(description="Count number of nucleotides for each sequence in a fasta file. Returns as genome file ready for input to software like bedtools. Prints to stdout.")

parser.add_argument("-f", "--fasta", type = str, required = True, 
            help = "Correctly formatted fasta.")

args = parser.parse_args()


def fasta_count(file_name):
    seq_dict = dict()
    with open(file_name) as fa:
        for ln in fa:
            line = ln.strip()
            if len(line) > 0:
                if line[0] == ">":
                    seq_name = line[1:]
                    if seq_name not in seq_dict:
                        seq_dict[seq_name] = 0
                    else:
                        raise ValueError("Fasta headers are not unique.")
                else:
                    seq_dict[seq_name] += len(ln.strip())
    return seq_dict

count_dict = fasta_count(args.fasta)

for header, count in count_dict.items():
    print(f"{header}\t{count}")

