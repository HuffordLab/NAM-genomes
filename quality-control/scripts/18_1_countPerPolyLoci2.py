"""
This script is designed to take a hapMap format SNP file and count percent
heterozygous loci based on the formula: (#heterozygotes / # all loci) excluding
missing data, average information over all loci, and report values per 
individual in an output file with 2 columns: indName  perHetLoci  
Created by David E. Hufnagel Jan 25, 2018

#This version is for standard hapmap format as opposed to 1 letter per diploid genotype
"""
import sys

inp = open(sys.argv[1])
out = open(sys.argv[2], "w")

        
def SaveIntoDict(key, val, dictx):
    if key not in dictx:
        dictx[key] = [val,]
    else:
        dictx[key].append(val)
        
def CalcPerHet(genos):
    hets = 0
    total = 0 #total excludes missing data
    for geno in genos:
        if geno != "NA":
            if geno == "het":
                hets += 1
            total += 1

    if total != 0:
        perHet = float(hets) / float(total)
        return perHet
    else:
        return "NA"


#Go through inp and gather SNP data in a dict of key: ind val: list of het or homo
#  as well as a list of individuals in the original order
title = inp.readline()
inds = title.strip().split("\t")[11:]
genoDict = {}
for line in inp:
    lineLst = line.strip().split("\t")
    genos = lineLst[11:]
    alleles = lineLst[1]    

    cnt = 0
    for geno in genos:
        ind = inds[cnt]
        if geno == "NA":
            SaveIntoDict(ind, "NA", genoDict)
        elif geno[0] != geno[1]:
            SaveIntoDict(ind, "het", genoDict)
        else:
            SaveIntoDict(ind, "homo", genoDict)
        cnt += 1

#Go through dict and count number of heterzygous loci and number of non-missing
#  data loci for each individual, storing the info in a new dict of 
#  key: ind val: percent heteryzygosity
perHetDict = {}
for key, val in genoDict.items():
    perHet = CalcPerHet(val)
    perHetDict[key] = perHet

#Use the percent heterozygosity dict and individuals list to create the output
for ind in inds:
    perHet = perHetDict[ind]
    newLine = "%s\t%s\n" % (ind, perHet)
    out.write(newLine)



inp.close()
out.close()

