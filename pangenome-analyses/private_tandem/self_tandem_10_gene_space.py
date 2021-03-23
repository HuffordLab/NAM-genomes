putative_self_tandem = []
with open("NAM_all_self_tandem.txt") as tandem_pair:
    for line in tandem_pair:
        gene_pair = line.strip().split(",") # Split each line.
        #print(gene_pair
        gene_space_diff = int(gene_pair[0][9:15])-int(gene_pair[1][9:15])
        if abs(gene_space_diff) < 11:
            tandem_pair = gene_pair[0] + "," + gene_pair[1] 
            putative_self_tandem.append(tandem_pair)
