"""
This script is designed to take a large hapmap format SNP file and subsample
based on a random probability of keeping each non indel locus.  indel loci
are all removed.  Also filters loci by missing data
Created by David E. Hufnagel on Jan 24, 2018
"""
import sys, random

inp = open(sys.argv[1])
out = open(sys.argv[2], "w")
thresh = float(sys.argv[3]) / 100.0       #The percent chance of keeping each marker
perMDthresh = float(sys.argv[4]) / 100.0  #The percent missing data allowed


#calculates percent missing data
def calcPerMD(genos):
    MDcnt = 0
    total = 0
    for geno in genos:
        if geno == "N":
            MDcnt += 1
        total += 1
        
    return (float(MDcnt) / float(total))
            

out.write(inp.readline())
for line in inp:
    lineLst = line.strip().split("\t")
    if lineLst[1] not in ["+/-","-/+"]: #filter out indels
        genos = lineLst[11:]
        perMD = calcPerMD(genos)
        
        #fliter for missing data
        if perMD <= perMDthresh:
            #random subsampling
            randNum = random.uniform(0,1)
            if randNum <= thresh:
                out.write(line)



inp.close()
out.close()
