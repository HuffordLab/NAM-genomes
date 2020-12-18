"""
This script is designed to remove problematic characters from fasta names
Created by David E. Hufnagel Jan 24, 2018
"""
import sys

fasta = open(sys.argv[1])
out = open(sys.argv[2], "w")

for line in fasta:
    if line.startswith(">"):
        name = line.strip().strip(">")
        newName = name.replace(" ","_").replace(":","_").replace(";","_").replace("[","_").replace("]","_").replace("(","_").replace(")","_")
        out.write(">" + newName + "\n")
    else:
        out.write(line)



fasta.close()
out.close()
