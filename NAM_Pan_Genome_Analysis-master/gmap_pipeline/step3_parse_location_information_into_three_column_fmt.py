gff_intersect = []
with open("/Users/yinjieqiu/Desktop/NAM_pan_genome_final/fill_in_info_include_scaf/new_gmap_gff_intersect/75_intersect_80_gmap_cutoff.gff") as gff3:
    for gff_line in gff3:
        # Split each line.
        # Split each line.
        gff_split = gff_line.strip().split()
        #print(gff_split)
        gmap_chr = gff_split[0]
        pan_gene_chr = gff_split[9]
        gmap_genome_name = gff_split[1]
        pan_genome_name = gff_split[10]
        gmap_stand_direction = gff_split[6]
        pan_gene_strand_direction = gff_split[15]
        #print(pan_genome_name)
        if gmap_chr == pan_gene_chr and gmap_genome_name == pan_genome_name and gmap_stand_direction == pan_gene_strand_direction:
            pan_gene_ID = gff_split[8].split(";")[0].split("=")[1]
            gmap_genome = gff_split[1]
            gmap_coordinate = gff_split[8].split(";")[1].split("=")[1]
            gene_model_associate = gff_split[17].split(";")[0].split("=")[1]
            gene_info = (pan_gene_ID + "_" + gmap_genome + "," + gene_model_associate + ";" + gmap_coordinate)
            #print(gene_info)
            gff_intersect.append(gene_info)
                    
                    
with open('/Users/yinjieqiu/Desktop/NAM_pan_genome_final/fill_in_info_include_scaf/new_gmap_gff_intersect/gff_intersect_cutoff_80.txt', 'w') as f:
    for item in gff_intersect:
        f.write("%s\n" % item)

        
# making dictionary for gff files to speed up the search 
# this will be the list made from above 
fill_in_info = {}
with open('/Users/yinjieqiu/Desktop/NAM_pan_genome_final/fill_in_info_include_scaf/new_gmap_gff_intersect/gff_intersect_cutoff_80.txt') as all_info:
    for line in all_info:
        # Split each line.
        split_fill_in_info = line.strip().split(",")
        #gene_info = [split_fill_in_info[1]]
        #this gene info could be changed into gene name only as below 
        gene_info = [split_fill_in_info[1].split(";")[0]]
        dict_keys_list = [split_fill_in_info[0]]
        d = dict(zip(dict_keys_list, gene_info))
        print(d)
        fill_in_info.update(d)
        
        
        
#update the list before input to R for merging 
# going through each lines in the fill in info file, if the first column is in the dictionary, then print first column and the value from the dictionary

final_list = []
with open("/Users/yinjieqiu/Desktop/NAM_pan_genome_final/fill_in_info_include_scaf/new_gmap_gff_intersect/gmap_info_for_fill_80.txt") as all_info:
    for line in all_info:
        line_split = line.strip().split()
        line_ID = line_split[0] + "_" + line_split[1]
        if line_ID in fill_in_info.keys(): 
            value = fill_in_info.get(line_ID)
            pan_gene_info = line_ID.split("_")[0] + "_" + line_ID.split("_")[1] + "\t" + line_ID.split("_")[2]
            final_info = pan_gene_info + "\t" + value
            print(final_info)
            final_list.append(final_info)
        else:
            #print(line_split)
            final_list.append(line)

with open('/Users/yinjieqiu/Desktop/NAM_pan_genome_final/fill_in_info_include_scaf/new_gmap_gff_intersect/final_list_for_merge_90.txt', 'w') as f:
    for item in final_list:
        f.write("%s\n" % item)
