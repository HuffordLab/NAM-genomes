"""
This script was designed to filter the HapMap2 SNP file in hapmap format such
that it contains only NAM founders
Created by David E. Hufnagel on Feb 19, 2018

Version two includes teosinte
"""
import sys
inp = open(sys.argv[1])
out = open(sys.argv[2], "w")

#Determine which indices are the NAM founders
title = inp.readline()
inds = title.strip().split("\t")[11:]
NAMlst = ["IL14H","MO18W","CML322","CML103","KI3","CML69","CML333","KY21",\
"TZI8","M162W","CML247","NC358","B73","CML277","TX303","B97","OH7B","P39",\
"CML52","OH43","MS71","HP301","M37W","CML228", "TIL11", "NC350", "KI11"]
cnt = 0
NAMindices = []
NAMinds = []
for ind in inds:
    ind = ind.split(":")[0]
    if ind in NAMlst:
        NAMindices.append(cnt)
        NAMinds.append(ind)
    cnt += 1   
     
#output title
newTitle = "\t".join(title.strip().split("\t")[:11]) + "\t" + "\t".join(NAMinds) + "\n"
out.write(newTitle)

#output each line for only NAM founders
for line in inp:
    lineLst = line.strip().split("\t")
    genos = lineLst[11:]
    cnt = 0
    newGenos = []
    for geno in genos:
        if cnt in NAMindices:
            newGenos.append(geno)
        cnt += 1

    newLine = "\t".join(lineLst[:11]) + "\t" + "\t".join(newGenos) + "\n"
    out.write(newLine)



inp.close()
out.close()
