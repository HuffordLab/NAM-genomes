import sys

sv = sys.argv[1]
with open(sv) as f:
    next(f)
    for line in f:
        ln = line.strip().split()
        #tra.chr7.29826131       T/A     1       25957
        #del.chr1.8840.12032
        svinfo_dot, genotype, chr2, pos2or3 =  ln[0:4]
        svinfo = svinfo_dot.split(".")
        if svinfo[0] == "tra":
            svtype, chr1, pos1 = svinfo
            if "scaf" not in chr1 and "scaf" not in chr2:
                print(f"{chr1}\t{int(pos1)-1}\t{pos1}\t{svinfo_dot}", flush = True)
                print(f"chr{chr2}\t{int(pos2)-1}\t{pos2}\t{svinfo_dot}", flush = True) 
        else:
            svtype, chr1, pos1, pos2 = svinfo
            if "scaf" not in chr1 and "scaf" not in chr2:
                print(f"{chr1}\t{int(pos1)-1}\t{pos2}\t{svinfo_dot}", flush = True)
        
