tandem_dict_for_sub ={}
with open("tandem_for_dict.txt") as tandem_pair:
    for line in tandem_pair:
        gene_pair = line.strip().split("\t") # Split each line.
        pan_gene_anchor = [gene_pair[0]]
        tandem_pair_info = [gene_pair[1]]
        d = dict(zip(pan_gene_anchor, tandem_pair_info))
        tandem_dict_for_sub.update(d)
        
final_list = []
with open("private_gene_matrix.csv") as private_gene_info:
    for line in private_gene_info:
        line_split = line.strip().split(",")
        if line_split[1] in tandem_dict_for_sub.keys():
            sub_line = line_split[2:35]
            for n, i in enumerate(sub_line):
                if i == line_split[1]:
                    sub_line[n] = tandem_dict_for_sub.get(line_split[1])
                    replace_line = line_split[0:2]+sub_line
                    final_list.append(replace_line)
        else: 
            final_list.append(line_split)

with open('tamden_add_private_test.txt', 'w') as f:
    for item in final_list:
        print_item = str(item).replace('[','').replace(']','')
        f.write("%s\n" % print_item)
