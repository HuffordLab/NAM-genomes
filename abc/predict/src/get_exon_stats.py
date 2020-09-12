#gff="data/ref/Zm-B73-REFERENCE-NAM-5.0_Zm00001e.1.gff"

#awk '$3 ~ /gene/' $gff | grep -v "scaf" | bedtools sort -i stdin | bedtools merge -i stdin | bedtools complement -g data/ref/Zm-B73-REFERENCE-NAM-5.0.gbed -i stdin


import pandas as pd
import numpy as np

qs = [0, 0.025, 0.025, 0.5, 0.75, 0.975, 1]

names = ["seqname","source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
gff = pd.read_csv("data/ref/Zm-B73-REFERENCE-NAM-5.0_Zm00001e.1.gff", sep = "\t", comment = "#", names = names)

chroms = [f"chr{i}" for i in range(1,11)]
gff = gff[gff.seqname.isin(chroms)]

#get distances between cds on same chromosome
gff_cds = gff[gff['feature'] == "CDS"]

same_seq =  np.array(gff_cds.iloc[1:]["seqname"]) == np.array(gff_cds.iloc[0:-1]["seqname"])
gaps = np.array(gff_cds.iloc[1:]["start"]) - np.array(gff_cds.iloc[0:-1]["end"])

print("gaps", gaps[same_seq])
print(np.quantile(gaps[same_seq], qs))


#get count of cds per gene
cds_count = np.array(gff_cds.groupby(["attribute"])["start"].count().reset_index()["start"])
print("cds count:", cds_count)
print("cds count:", np.quantile(cds_count, qs))


#get distribution of cds sizes
print("cds sizes", np.array(gff_cds["end"] - gff_cds["start"]))
print("cds sizes", np.quantile(np.array(gff_cds["end"] - gff_cds["start"]), qs))

#get distances between gene features on the same chromosome
gff_gene = gff[gff['feature'] == "gene"]
same_seq =   np.array(gff_gene.iloc[1:]["seqname"]) == np.array(gff_gene.iloc[0:-1]["seqname"])
gap_gene = np.array(gff_gene.iloc[1:]["start"]) - np.array(gff_gene.iloc[0:-1]["end"])
print("between genes", gap_gene[same_seq])
print(np.quantile(gap_gene[same_seq], qs))

