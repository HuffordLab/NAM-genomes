"""
This script is designed to take a hapmap format SNP file as input and create
an output with the same format without monomorphic loci

Created by David E. Hufnagel on Feb 13, 2018
"""
import sys

inp = open(sys.argv[1])
out = open(sys.argv[2], "w")


def IsMonoLoci(snps):
    isMono = True
    firstLoc = ""
    for snp in snps:
        if snp == "N": #assumed to be the only missing data condition
            continue
        elif snp in ["R","Y","S","W","K","M"]: #a heterozygotic locus proves its not monomorphic
            isMono = False
            break
        elif firstLoc == "": #first homozygote
            firstLoc = snp
        else:
            if snp != firstLoc:
                isMono = False
                break

    return(isMono)

def SaveIntoDict(key, val, dictx):
    if not key in dictx:
        dictx[key] = [val,]
    else:
        dictx[key].append(val)

#takes a list of lists and inverts the matrix
def InvertMatrix(matrix):
    newMatrixDict = {}
    for col in matrix:
        cnt = 0
        for row in col:
            SaveIntoDict(cnt, row, newMatrixDict)
            cnt += 1

    cnt = 0
    newMatrix = []
    while cnt in newMatrixDict:
        data = newMatrixDict[cnt]
        newMatrix.append(data)
        cnt += 1

    return(newMatrix)


out.write(inp.readline())
mml = 0 #count of monomorphic loci
indels = 0
for line in inp:
    lineLst = line.strip().split()
    alleles = lineLst[1]
    
    #skip indels
    if not alleles in ["+/-","-/+"]:
        genos = lineLst[11:]
        isMono = IsMonoLoci(genos)
        if isMono:
            mml += 1
        else:
            out.write(line)
    else:
        indels += 1

print "%s indels removed" % (indels)  
print "%s monomorphic loci removed" % (mml)        


inp.close()
out.close()