"""
This script was created to expand a condensed hapmap format SNP file with 1
character representing one locus for a diploid individual to the standard 
format with 2 characters per locus per individual.
Created by David E. Hufnagel on Feb 13, 2018
"""
import sys

inp = open(sys.argv[1])
out = open(sys.argv[2], "w")


def ExpandGeno(old):
    if old == "A":
        new = "AA"
    elif old == "T":
        new = "TT"
    elif old == "C":
        new = "CC"
    elif old == "G":
        new = "GG"
    elif old == "N":
        new = "NA"
    elif old == "R":
        new = "AG"
    elif old == "Y":
        new = "CT"
    elif old == "S":
        new = "GC"
    elif old == "W":
        new = "AT"
    elif old == "K":
        new = "GT"
    elif old == "M":
        new = "AC"
    else:
        print "ERROR"
        print old
        sys.exit()
        
    return(new)


out.write(inp.readline())
for line in inp:
    lineLst = line.strip().split()
    alleles = lineLst[1]
    
    #skip indels
    if not alleles in ["+/-","-/+"]:
        genos = lineLst[11:]
        newGenos = []
        for geno in genos:
            newGeno = ExpandGeno(geno)
            newGenos.append(newGeno)

        allData = lineLst[:11]
        allData.extend(newGenos)
        newLine = "\t".join(allData) + "\n"
        out.write(newLine)
        
        

inp.close()
out.close()