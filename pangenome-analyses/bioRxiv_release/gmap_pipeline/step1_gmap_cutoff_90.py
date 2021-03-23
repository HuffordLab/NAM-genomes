# making dictionary for gff files to speed up the search 
master_gff = {}
with open("/home/hirschc1/qiuxx221/nam_pan_genome/NAM_annotation/canonical_gff/modified_fmt_pan_gene.gff") as gff3:
    for gff_line in gff3:
        # Split each line.
        gff_line = gff_line.strip().split()
        gene_info = [gff_line[0]+ "_" + gff_line[1] + "_" + gff_line[3]+ "_" + gff_line[4] + "_" + gff_line[6]]
        gene_ID = gff_line[8].split(";")[0].split("=")[1]
        dict_keys_list = [gene_ID]
        d = dict(zip(dict_keys_list, gene_info))
        master_gff.update(d)
        
# making gff dict for all NAM genomes     
B73_dict = {}
with open("/home/hirschc1/qiuxx221/nam_pan_genome/nucmer_1000_subgenome/pan_to_B73_1_path.gff3") as gff3:
    for gff_line in gff3:
        # Split each line.
        gff_line = gff_line.strip().split()
        sequence_coverage = float(gff_line[8].split(";")[3].split("=")[1])
        sequence_identity = float(gff_line[8].split(";")[4].split("=")[1]) 
        if sequence_coverage > 90 and sequence_identity >90: 
            gene_info = [gff_line[0]+ "_" + gff_line[1] + "_" + gff_line[3]+ "_" + gff_line[4] + "_" + gff_line[6]]
            gene_ID = gff_line[8].split(";")[1].split("=")[1]
            dict_keys_list = [gene_ID]
            d = dict(zip(dict_keys_list, gene_info))
            B73_dict.update(d)
        
B97_dict = {}
with open("/home/hirschc1/qiuxx221/nam_pan_genome/nucmer_1000_subgenome/pan_to_B97_1_path.gff3") as gff3:
    for gff_line in gff3:
        # Split each line.
        gff_line = gff_line.strip().split()
        sequence_coverage = float(gff_line[8].split(";")[3].split("=")[1])
        sequence_identity = float(gff_line[8].split(";")[4].split("=")[1]) 
        if sequence_coverage > 90 and sequence_identity >90: 
            gene_info = [gff_line[0]+ "_" + gff_line[1] + "_" + gff_line[3]+ "_" + gff_line[4] + "_" + gff_line[6]]
            gene_ID = gff_line[8].split(";")[1].split("=")[1]
            dict_keys_list = [gene_ID]
            d = dict(zip(dict_keys_list, gene_info))
            B97_dict.update(d)
            
CML52_dict = {}
with open("/home/hirschc1/qiuxx221/nam_pan_genome/nucmer_1000_subgenome/pan_to_CML52_1_path.gff3") as gff3:
    for gff_line in gff3:
        # Split each line.
        gff_line = gff_line.strip().split()
        sequence_coverage = float(gff_line[8].split(";")[3].split("=")[1])
        sequence_identity = float(gff_line[8].split(";")[4].split("=")[1]) 
        if sequence_coverage > 90 and sequence_identity >90: 
            gene_info = [gff_line[0]+ "_" + gff_line[1] + "_" + gff_line[3]+ "_" + gff_line[4] + "_" + gff_line[6]]
            gene_ID = gff_line[8].split(";")[1].split("=")[1]
            dict_keys_list = [gene_ID]
            d = dict(zip(dict_keys_list, gene_info))
            CML52_dict.update(d)
            
CML69_dict = {}
with open("/home/hirschc1/qiuxx221/nam_pan_genome/nucmer_1000_subgenome/pan_to_CML69_1_path.gff3") as gff3:
    for gff_line in gff3:
        # Split each line.
        gff_line = gff_line.strip().split()
        sequence_coverage = float(gff_line[8].split(";")[3].split("=")[1])
        sequence_identity = float(gff_line[8].split(";")[4].split("=")[1]) 
        if sequence_coverage > 90 and sequence_identity > 90: 
            gene_info = [gff_line[0]+ "_" + gff_line[1] + "_" + gff_line[3]+ "_" + gff_line[4] + "_" + gff_line[6]]
            gene_ID = gff_line[8].split(";")[1].split("=")[1]
            dict_keys_list = [gene_ID]
            d = dict(zip(dict_keys_list, gene_info))
            CML69_dict.update(d)

CML103_dict = {}
with open("/home/hirschc1/qiuxx221/nam_pan_genome/nucmer_1000_subgenome/pan_to_CML103_1_path.gff3") as gff3:
    for gff_line in gff3:
        # Split each line.
        gff_line = gff_line.strip().split()
        sequence_coverage = float(gff_line[8].split(";")[3].split("=")[1])
        sequence_identity = float(gff_line[8].split(";")[4].split("=")[1]) 
        if sequence_coverage > 90 and sequence_identity > 90: 
            gene_info = [gff_line[0]+ "_" + gff_line[1] + "_" + gff_line[3]+ "_" + gff_line[4] + "_" + gff_line[6]]
            gene_ID = gff_line[8].split(";")[1].split("=")[1]
            dict_keys_list = [gene_ID]
            d = dict(zip(dict_keys_list, gene_info))
            CML103_dict.update(d)
        
CML228_dict = {}
with open("/home/hirschc1/qiuxx221/nam_pan_genome/nucmer_1000_subgenome/pan_to_CML228_1_path.gff3") as gff3:
    for gff_line in gff3:
        # Split each line.
        gff_line = gff_line.strip().split()
        sequence_coverage = float(gff_line[8].split(";")[3].split("=")[1])
        sequence_identity = float(gff_line[8].split(";")[4].split("=")[1]) 
        if sequence_coverage > 90 and sequence_identity > 90: 
            gene_info = [gff_line[0]+ "_" + gff_line[1] + "_" + gff_line[3]+ "_" + gff_line[4] + "_" + gff_line[6]]
            gene_ID = gff_line[8].split(";")[1].split("=")[1]
            dict_keys_list = [gene_ID]
            d = dict(zip(dict_keys_list, gene_info))
            CML228_dict.update(d)

CML247_dict = {}
with open("/home/hirschc1/qiuxx221/nam_pan_genome/nucmer_1000_subgenome/pan_to_CML247_1_path.gff3") as gff3:
    for gff_line in gff3:
        # Split each line.
        gff_line = gff_line.strip().split()
        sequence_coverage = float(gff_line[8].split(";")[3].split("=")[1])
        sequence_identity = float(gff_line[8].split(";")[4].split("=")[1]) 
        if sequence_coverage > 90 and sequence_identity > 90: 
            gene_info = [gff_line[0]+ "_" + gff_line[1] + "_" + gff_line[3]+ "_" + gff_line[4] + "_" + gff_line[6]]
            gene_ID = gff_line[8].split(";")[1].split("=")[1]
            dict_keys_list = [gene_ID]
            d = dict(zip(dict_keys_list, gene_info))
            CML247_dict.update(d)
        
CML277_dict = {}
with open("/home/hirschc1/qiuxx221/nam_pan_genome/nucmer_1000_subgenome/pan_to_CML277_1_path.gff3") as gff3:
    for gff_line in gff3:
        # Split each line.
        gff_line = gff_line.strip().split()
        sequence_coverage = float(gff_line[8].split(";")[3].split("=")[1])
        sequence_identity = float(gff_line[8].split(";")[4].split("=")[1]) 
        if sequence_coverage > 90 and sequence_identity > 90: 
            gene_info = [gff_line[0]+ "_" + gff_line[1] + "_" + gff_line[3]+ "_" + gff_line[4] + "_" + gff_line[6]]
            gene_ID = gff_line[8].split(";")[1].split("=")[1]
            dict_keys_list = [gene_ID]
            d = dict(zip(dict_keys_list, gene_info))
            CML277_dict.update(d)

CML322_dict = {}
with open("/home/hirschc1/qiuxx221/nam_pan_genome/nucmer_1000_subgenome/pan_to_CML322_1_path.gff3") as gff3:
    for gff_line in gff3:
        # Split each line.
        gff_line = gff_line.strip().split()
        sequence_coverage = float(gff_line[8].split(";")[3].split("=")[1])
        sequence_identity = float(gff_line[8].split(";")[4].split("=")[1]) 
        if sequence_coverage > 90 and sequence_identity > 90: 
            gene_info = [gff_line[0]+ "_" + gff_line[1] + "_" + gff_line[3]+ "_" + gff_line[4] + "_" + gff_line[6]]
            gene_ID = gff_line[8].split(";")[1].split("=")[1]
            dict_keys_list = [gene_ID]
            d = dict(zip(dict_keys_list, gene_info))
            CML322_dict.update(d)
        

CML333_dict = {}
with open("/home/hirschc1/qiuxx221/nam_pan_genome/nucmer_1000_subgenome/pan_to_CML333_1_path.gff3") as gff3:
    for gff_line in gff3:
        # Split each line.
        gff_line = gff_line.strip().split()
        sequence_coverage = float(gff_line[8].split(";")[3].split("=")[1])
        sequence_identity = float(gff_line[8].split(";")[4].split("=")[1]) 
        if sequence_coverage > 90 and sequence_identity > 90: 
            gene_info = [gff_line[0]+ "_" + gff_line[1] + "_" + gff_line[3]+ "_" + gff_line[4] + "_" + gff_line[6]]
            gene_ID = gff_line[8].split(";")[1].split("=")[1]
            dict_keys_list = [gene_ID]
            d = dict(zip(dict_keys_list, gene_info))
            CML333_dict.update(d)
        
HP301_dict = {}
with open("/home/hirschc1/qiuxx221/nam_pan_genome/nucmer_1000_subgenome/pan_to_HP301_1_path.gff3") as gff3:
    for gff_line in gff3:
        # Split each line.
        gff_line = gff_line.strip().split()
        sequence_coverage = float(gff_line[8].split(";")[3].split("=")[1])
        sequence_identity = float(gff_line[8].split(";")[4].split("=")[1]) 
        if sequence_coverage > 90 and sequence_identity > 90: 
            gene_info = [gff_line[0]+ "_" + gff_line[1] + "_" + gff_line[3]+ "_" + gff_line[4] + "_" + gff_line[6]]
            gene_ID = gff_line[8].split(";")[1].split("=")[1]
            dict_keys_list = [gene_ID]
            d = dict(zip(dict_keys_list, gene_info))
            HP301_dict.update(d)   
        
Il14H_dict = {}
with open("/home/hirschc1/qiuxx221/nam_pan_genome/nucmer_1000_subgenome/pan_to_Il14H_1_path.gff3") as gff3:
    for gff_line in gff3:
        # Split each line.
        gff_line = gff_line.strip().split()
        sequence_coverage = float(gff_line[8].split(";")[3].split("=")[1])
        sequence_identity = float(gff_line[8].split(";")[4].split("=")[1]) 
        if sequence_coverage > 90 and sequence_identity > 90: 
            gene_info = [gff_line[0]+ "_" + gff_line[1] + "_" + gff_line[3]+ "_" + gff_line[4] + "_" + gff_line[6]]
            gene_ID = gff_line[8].split(";")[1].split("=")[1]
            dict_keys_list = [gene_ID]
            d = dict(zip(dict_keys_list, gene_info))
            Il14H_dict.update(d)   

Ki3_dict = {}
with open("/home/hirschc1/qiuxx221/nam_pan_genome/nucmer_1000_subgenome/pan_to_Ki3_1_path.gff3") as gff3:
    for gff_line in gff3:
        # Split each line.
        gff_line = gff_line.strip().split()
        sequence_coverage = float(gff_line[8].split(";")[3].split("=")[1])
        sequence_identity = float(gff_line[8].split(";")[4].split("=")[1]) 
        if sequence_coverage > 90 and sequence_identity > 90: 
            gene_info = [gff_line[0]+ "_" + gff_line[1] + "_" + gff_line[3]+ "_" + gff_line[4] + "_" + gff_line[6]]
            gene_ID = gff_line[8].split(";")[1].split("=")[1]
            dict_keys_list = [gene_ID]
            d = dict(zip(dict_keys_list, gene_info))
            Ki3_dict.update(d)   
        
Ki11_dict = {}
with open("/home/hirschc1/qiuxx221/nam_pan_genome/nucmer_1000_subgenome/pan_to_Ki11_1_path.gff3") as gff3:
    for gff_line in gff3:
        # Split each line.
        gff_line = gff_line.strip().split()
        sequence_coverage = float(gff_line[8].split(";")[3].split("=")[1])
        sequence_identity = float(gff_line[8].split(";")[4].split("=")[1]) 
        if sequence_coverage > 90 and sequence_identity > 90: 
            gene_info = [gff_line[0]+ "_" + gff_line[1] + "_" + gff_line[3]+ "_" + gff_line[4] + "_" + gff_line[6]]
            gene_ID = gff_line[8].split(";")[1].split("=")[1]
            dict_keys_list = [gene_ID]
            d = dict(zip(dict_keys_list, gene_info))
            Ki11_dict.update(d)         
        
M162W_dict = {}
with open("/home/hirschc1/qiuxx221/nam_pan_genome/nucmer_1000_subgenome/pan_to_M162W_1_path.gff3") as gff3:
    for gff_line in gff3:
        # Split each line.
        # Split each line.
        gff_line = gff_line.strip().split()
        sequence_coverage = float(gff_line[8].split(";")[3].split("=")[1])
        sequence_identity = float(gff_line[8].split(";")[4].split("=")[1]) 
        if sequence_coverage > 90 and sequence_identity > 90: 
            gene_info = [gff_line[0]+ "_" + gff_line[1] + "_" + gff_line[3]+ "_" + gff_line[4] + "_" + gff_line[6]]
            gene_ID = gff_line[8].split(";")[1].split("=")[1]
            dict_keys_list = [gene_ID]
            d = dict(zip(dict_keys_list, gene_info))
            M162W_dict.update(d)     
        
M37W_dict = {}
with open("/home/hirschc1/qiuxx221/nam_pan_genome/nucmer_1000_subgenome/pan_to_M37W_1_path.gff3") as gff3:
    for gff_line in gff3:
        # Split each line.
        gff_line = gff_line.strip().split()
        sequence_coverage = float(gff_line[8].split(";")[3].split("=")[1])
        sequence_identity = float(gff_line[8].split(";")[4].split("=")[1]) 
        if sequence_coverage > 90 and sequence_identity > 90: 
            gene_info = [gff_line[0]+ "_" + gff_line[1] + "_" + gff_line[3]+ "_" + gff_line[4] + "_" + gff_line[6]]
            gene_ID = gff_line[8].split(";")[1].split("=")[1]
            dict_keys_list = [gene_ID]
            d = dict(zip(dict_keys_list, gene_info))
            M37W_dict.update(d) 
        
Mo18W_dict = {}
with open("/home/hirschc1/qiuxx221/nam_pan_genome/nucmer_1000_subgenome/pan_to_Mo18W_1_path.gff3") as gff3:
    for gff_line in gff3:
        # Split each line.
        gff_line = gff_line.strip().split()
        sequence_coverage = float(gff_line[8].split(";")[3].split("=")[1])
        sequence_identity = float(gff_line[8].split(";")[4].split("=")[1]) 
        if sequence_coverage > 90 and sequence_identity > 90: 
            gene_info = [gff_line[0]+ "_" + gff_line[1] + "_" + gff_line[3]+ "_" + gff_line[4] + "_" + gff_line[6]]
            gene_ID = gff_line[8].split(";")[1].split("=")[1]
            dict_keys_list = [gene_ID]
            d = dict(zip(dict_keys_list, gene_info))
            Mo18W_dict.update(d)         
        
NC350_dict = {}
with open("/home/hirschc1/qiuxx221/nam_pan_genome/nucmer_1000_subgenome/pan_to_NC350_1_path.gff3") as gff3:
    for gff_line in gff3:
        # Split each line.
        gff_line = gff_line.strip().split()
        sequence_coverage = float(gff_line[8].split(";")[3].split("=")[1])
        sequence_identity = float(gff_line[8].split(";")[4].split("=")[1]) 
        if sequence_coverage > 90 and sequence_identity > 90: 
            gene_info = [gff_line[0]+ "_" + gff_line[1] + "_" + gff_line[3]+ "_" + gff_line[4] + "_" + gff_line[6]]
            gene_ID = gff_line[8].split(";")[1].split("=")[1]
            dict_keys_list = [gene_ID]
            d = dict(zip(dict_keys_list, gene_info))
            NC350_dict.update(d)           
        
NC358_dict = {}
with open("/home/hirschc1/qiuxx221/nam_pan_genome/nucmer_1000_subgenome/pan_to_NC358_1_path.gff3") as gff3:
    for gff_line in gff3:
        # Split each line.
        gff_line = gff_line.strip().split()
        sequence_coverage = float(gff_line[8].split(";")[3].split("=")[1])
        sequence_identity = float(gff_line[8].split(";")[4].split("=")[1]) 
        if sequence_coverage > 90 and sequence_identity > 90: 
            gene_info = [gff_line[0]+ "_" + gff_line[1] + "_" + gff_line[3]+ "_" + gff_line[4] + "_" + gff_line[6]]
            gene_ID = gff_line[8].split(";")[1].split("=")[1]
            dict_keys_list = [gene_ID]
            d = dict(zip(dict_keys_list, gene_info))
            NC358_dict.update(d)            

Oh7B_dict = {}
with open("/home/hirschc1/qiuxx221/nam_pan_genome/nucmer_1000_subgenome/pan_to_Oh7B_1_path.gff3") as gff3:
    for gff_line in gff3:
        # Split each line.
        gff_line = gff_line.strip().split()
        sequence_coverage = float(gff_line[8].split(";")[3].split("=")[1])
        sequence_identity = float(gff_line[8].split(";")[4].split("=")[1]) 
        if sequence_coverage > 90 and sequence_identity > 90: 
            gene_info = [gff_line[0]+ "_" + gff_line[1] + "_" + gff_line[3]+ "_" + gff_line[4] + "_" + gff_line[6]]
            gene_ID = gff_line[8].split(";")[1].split("=")[1]
            dict_keys_list = [gene_ID]
            d = dict(zip(dict_keys_list, gene_info))
            Oh7B_dict.update(d)          
        
Oh43_dict = {}
with open("/home/hirschc1/qiuxx221/nam_pan_genome/nucmer_1000_subgenome/pan_to_Oh43_1_path.gff3") as gff3:
    for gff_line in gff3:
        # Split each line.
        gff_line = gff_line.strip().split()
        sequence_coverage = float(gff_line[8].split(";")[3].split("=")[1])
        sequence_identity = float(gff_line[8].split(";")[4].split("=")[1]) 
        if sequence_coverage > 90 and sequence_identity > 90: 
            gene_info = [gff_line[0]+ "_" + gff_line[1] + "_" + gff_line[3]+ "_" + gff_line[4] + "_" + gff_line[6]]
            gene_ID = gff_line[8].split(";")[1].split("=")[1]
            dict_keys_list = [gene_ID]
            d = dict(zip(dict_keys_list, gene_info))
            Oh43_dict.update(d)      
        
P39_dict = {}
with open("/home/hirschc1/qiuxx221/nam_pan_genome/nucmer_1000_subgenome/pan_to_P39_1_path.gff3") as gff3:
    for gff_line in gff3:
        # Split each line.
        # Split each line.
        gff_line = gff_line.strip().split()
        sequence_coverage = float(gff_line[8].split(";")[3].split("=")[1])
        sequence_identity = float(gff_line[8].split(";")[4].split("=")[1]) 
        if sequence_coverage > 90 and sequence_identity > 90: 
            gene_info = [gff_line[0]+ "_" + gff_line[1] + "_" + gff_line[3]+ "_" + gff_line[4] + "_" + gff_line[6]]
            gene_ID = gff_line[8].split(";")[1].split("=")[1]
            dict_keys_list = [gene_ID]
            d = dict(zip(dict_keys_list, gene_info))
            P39_dict.update(d) 
        
Ky21_dict = {}
with open("/home/hirschc1/qiuxx221/nam_pan_genome/nucmer_1000_subgenome/pan_to_Ky21_1_path.gff3") as gff3:
    for gff_line in gff3:
        # Split each line.
        gff_line = gff_line.strip().split()
        sequence_coverage = float(gff_line[8].split(";")[3].split("=")[1])
        sequence_identity = float(gff_line[8].split(";")[4].split("=")[1]) 
        if sequence_coverage > 90 and sequence_identity > 90: 
            gene_info = [gff_line[0]+ "_" + gff_line[1] + "_" + gff_line[3]+ "_" + gff_line[4] + "_" + gff_line[6]]
            gene_ID = gff_line[8].split(";")[1].split("=")[1]
            dict_keys_list = [gene_ID]
            d = dict(zip(dict_keys_list, gene_info))
            Ky21_dict.update(d)         

Ms71_dict = {}
with open("/home/hirschc1/qiuxx221/nam_pan_genome/nucmer_1000_subgenome/pan_to_Ms71_1_path.gff3") as gff3:
    for gff_line in gff3:
        gff_line = gff_line.strip().split()
        sequence_coverage = float(gff_line[8].split(";")[3].split("=")[1])
        sequence_identity = float(gff_line[8].split(";")[4].split("=")[1]) 
        if sequence_coverage > 90 and sequence_identity > 90: 
            gene_info = [gff_line[0]+ "_" + gff_line[1] + "_" + gff_line[3]+ "_" + gff_line[4] + "_" + gff_line[6]]
            gene_ID = gff_line[8].split(";")[1].split("=")[1]
            dict_keys_list = [gene_ID]
            d = dict(zip(dict_keys_list, gene_info))
            Ms71_dict.update(d) 
        
Tx303_dict = {}
with open("/home/hirschc1/qiuxx221/nam_pan_genome/nucmer_1000_subgenome/pan_to_Tx303_1_path.gff3") as gff3:
    for gff_line in gff3:
        # Split each line.
        gff_line = gff_line.strip().split()
        sequence_coverage = float(gff_line[8].split(";")[3].split("=")[1])
        sequence_identity = float(gff_line[8].split(";")[4].split("=")[1]) 
        if sequence_coverage > 90 and sequence_identity > 90: 
            gene_info = [gff_line[0]+ "_" + gff_line[1] + "_" + gff_line[3]+ "_" + gff_line[4] + "_" + gff_line[6]]
            gene_ID = gff_line[8].split(";")[1].split("=")[1]
            dict_keys_list = [gene_ID]
            d = dict(zip(dict_keys_list, gene_info))
            Tx303_dict.update(d)         
        
Tzi8_dict = {}
with open("/home/hirschc1/qiuxx221/nam_pan_genome/nucmer_1000_subgenome/pan_to_Tzi8_1_path.gff3") as gff3:
    for gff_line in gff3:
        gff_line = gff_line.strip().split()
        sequence_coverage = float(gff_line[8].split(";")[3].split("=")[1])
        sequence_identity = float(gff_line[8].split(";")[4].split("=")[1]) 
        if sequence_coverage > 90 and sequence_identity > 90: 
            gene_info = [gff_line[0]+ "_" + gff_line[1] + "_" + gff_line[3]+ "_" + gff_line[4] + "_" + gff_line[6]]
            gene_ID = gff_line[8].split(";")[1].split("=")[1]
            dict_keys_list = [gene_ID]
            d = dict(zip(dict_keys_list, gene_info))
            Tzi8_dict.update(d)          
            
# make all dict into a dict for future look up 
dict_name_value = [B73_dict, B97_dict, CML52_dict, CML69_dict, CML103_dict ,CML228_dict, CML247_dict, CML277_dict, CML322_dict , CML333_dict, HP301_dict, Il14H_dict, Ki3_dict, Ki11_dict, M162W_dict, M37W_dict, Mo18W_dict,NC350_dict,NC358_dict, Oh7B_dict, Oh43_dict, P39_dict, Ky21_dict, Ms71_dict, Tx303_dict,Tzi8_dict]
dict_name_list = ['B73',  'B97', 'CML52', 'CML69', 'CML103' ,'CML228', 'CML247', 'CML277', 'CML322' , 'CML333', 'HP301', 'Il14H', 'Ki3', 'Ki11', 'M162W', 'M37W', 'Mo18W', 'NC350', 'NC358', 'Oh7B', 'Oh43', 'P39', 'Ky21', 'Ms71', 'Tx303', 'Tzi8']
universal_dict = dict(zip(dict_name_list,dict_name_value))

import pandas as pd
# getting a list of pan gene that has NA in the matrix 
# information include pan gene ID, genome it is coming from, and genome where it was missing
infile = '/home/hirschc1/qiuxx221/nam_pan_genome/gmap_pipeline/pan26_all.collapsed.csv'
pan_matrix = open(infile, "r")

# get header info
header = pan_matrix.readline()
header = header.strip()
header = header.split(",")
# key for pan genome name 

id_list = []
for line in pan_matrix:
    pan_gene = line.strip()
    pan_gene = line.split(',')
    for idx, gene_id in enumerate(pan_gene): 
        #print(gene_id)
        if gene_id == "NA":
            gene_info = [pan_gene[1].split(';')[0], pan_gene[1][5:9], header[idx]]
            pan_gene_name= ["B73","Ki11","P39","NC350","B97","CML103","CML228","CML247","CML277","CML322","CML333","CML52","CML69","HP301","Il14H","Ki3","M162W","M37W","Mo18W","NC358","Oh7B","Oh43","Ky21","Tx303","Tzi8","Ms71"]
            pan_gene_ID = ["01eb","30ab",'40ab','36ab','18ab','21ab','22ab','23ab','24ab','25ab','26ab','19ab','20ab','27ab','28ab','29ab','33ab','32ab','34ab','37ab','38ab','39ab','31ab','41ab','42ab','35ab']
            pan_gene_key = dict(zip(pan_gene_ID,pan_gene_name))
            pan_genome_name = pan_gene_key[gene_info[1]]
            info_for_look_up = [gene_info[0],pan_genome_name,gene_info[2]]
            id_list.append(info_for_look_up)
        if gene_id != "NA":
            continue           

coord_list=[]
for line in id_list:
    gmap_gff_dict= universal_dict[line[2]]
    if line[0] in gmap_gff_dict:
        gmap_coordinate = gmap_gff_dict[line[0]]
        master_coordinate = master_gff[line[0]]
        output_info = line[0] +" "+ gmap_coordinate + " "+ master_coordinate
        #print(output_info)
        coord_list.append(output_info)
        
syn_gene_info_for_fill = []        
gmap_syn_gene_check = []
gene_to_fill_gff = []
syn_gene_info_for_matrix = []
for line in coord_list:
    query_meta=[]
    query_chr = line.split(" ")[1].split("_")[0]
    query_genome_name = line.split(" ")[1].split("_")[1]
    query_start_position = line.split(" ")[1].split("_")[2]
    query_end_position = line.split(" ")[1].split("_")[3]
    strand_direction = line.split(" ")[1].split("_")[4]
    query_info = [query_chr, query_start_position, query_end_position,strand_direction]
    query_meta.append(list(query_info))
    #print(line)
    
    #pan_reference
    pan_meta=[]
    pan_gene_name = line.split(" ")[0]
    pan_chr = line.split(" ")[2].split("_")[0]
    pan_genome_name = line.split(" ")[2].split("_")[1]
    pan_start_position = line.split(" ")[2].split("_")[2]
    pan_end_position = line.split(" ")[2].split("_")[3]
    pan_info = [pan_gene_name, pan_chr, pan_start_position,pan_end_position]
    pan_meta.append(pan_info)
    #print(pan_meta)
    if query_chr == pan_chr: 
        nucmer_file_name_for_look_up= "/home/hirschc1/qiuxx221/nucmer_1000_filter/" + pan_genome_name + "_" + query_genome_name   + "_c1000.fil.coords"
        info_line = [line, nucmer_file_name_for_look_up]
        gmap_syn_gene_check.append(info_line)
        #print(gmap_syn_gene_check)
        df=pd.read_csv(nucmer_file_name_for_look_up,index_col=None,low_memory=False,sep="\t")
        try:
            Dframe_index = df.index[(df["TAGS"] == pan_chr) & (df['S1'] <= int(pan_start_position))].tolist()[-1]
            Syn = df.iloc[(Dframe_index)].tolist()
            Two_sided = True 
        except IndexError as e:
                print("type error: ",
                      str(e),
                      "for gene",
                      pan_gene_name,
                      "at begining of chr ",
                      "not in syntenic region")
                Two_sided = False
                continue
                
        Q_Syn_region_start = int(Syn[0])
        Q_Syn_region_stop = int(Syn[1])
        Q_Syn_region_chrom = Syn[11]
        Syn_region_start = int(Syn[2])
        Syn_region_stop = int(Syn[3])
        Syn_region_chrom = Syn[12]

        Syn_Meta = [Q_Syn_region_chrom,
                        Q_Syn_region_start,
                        Q_Syn_region_stop,
                        Syn_region_chrom,
                        Syn_region_start,
                        Syn_region_stop,]
        #print(Syn_Meta)
       
        for query_line in query_meta: 
                chrom = query_line[0]
                start = int(query_line[1])
                stop = int(query_line[2])
                #print(query_line)
                if chrom == Syn_region_chrom:
                    if ((Syn_region_start)-500000) <= start <= \
                       ((Syn_region_stop) + 500000) and \
                       ((Syn_region_start)-500000) <= stop <= \
                       ((Syn_region_stop) + 500000):
                        syn_gene_info = ("pan_gene_ID=" + pan_gene_name + ";" + "gmap_ID=" +query_line[0]+ ":"+ query_line[1] + "-"+ query_line[2] + ";" + "Pan_gene_origin=" + pan_genome_name )
                        syn_gene_matrix = (pan_gene_name + "\t" + query_genome_name + "\t" + "gmap_ID=" +query_line[0]+ ":"+ query_line[1] + "-"+ query_line[2]  )
                        syn_gene_info_for_matrix.append(syn_gene_matrix)
                        #print(syn_gene_info)
                        syn_gene_gff =  (chrom + '\t' + query_genome_name + '\t' + "mRNA"  + '\t' + str(start) + '\t' + str(stop) + '\t' + "." + '\t' + str(query_line[3]) + '\t' + "." + '\t' + str(syn_gene_info))
                        gene_to_fill_gff.append(syn_gene_gff)
                        print(syn_gene_gff)
                        
with open('gmap_cutoff_90.gff', 'w') as f:
    for item in gene_to_fill_gff:
        f.write("%s\n" % item)
with open('gmap_info_for_fill_90.txt', 'w') as f:
    for item in syn_gene_info_for_matrix:
        f.write("%s\n" % item)            
