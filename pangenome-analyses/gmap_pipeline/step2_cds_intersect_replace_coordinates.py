# master list of CDS gffusing gene annotation file (this includes both canonical and other transcripts)
cat *.gff | grep CDS > master_CDS.gff 
/home/hirschc1/qiuxx221/nam_pan_genome/NAM_annotation/gff/master_CDS.gff 


# coverage and identity cutoff ID information and get ID for  
cut -f 1 gmap_info_for_fill_90.txt > cov_id_90_ID.txt
/home/hirschc1/qiuxx221/nam_pan_genome/gmap_pipeline/cov_id_90_ID.txt

# master gmap CDS gff
cat pan_to_* | grep CDS > master_gmap_CDS.gff

# get cds that meet this standard
grep -Fwf /home/hirschc1/qiuxx221/nam_pan_genome/gmap_pipeline/cov_id_90_ID.txt master_gmap_CDS.gff > 9090_cutoff_CDS.gff

# 2851973 lines 

# filter those gmap coordinates with cds length more than 200 bp 
awk '($5-$4> 200)' 9090_cutoff_CDS.gff > 9090_CDS_200.gff
# 886902 lines 

# intersect those gmap coordiantes with gff annotation cds gff files, that requires reciprocal of 90% coverage at intersection 
intersectBed -a 9090_CDS_200.gff -b /home/hirschc1/qiuxx221/nam_pan_genome/NAM_annotation/gff/master_CDS.gff -r -f 0.90 -wa -wb > CDS_909090_gmap_cutoff.gff
# 240367 lines have intersect. 

# parse out intersection gff in python

gff_intersect = []
with open("/Users/yinjieqiu/Desktop/pan_genome_nov 2/CDS_909090_gmap_cutoff.gff") as gff3:
    for gff_line in gff3:
        # Split each line.
        # Split each line.
        gff_split = gff_line.strip().split()
        #print(gff_split)
        gmap_chr = gff_split[0]
        pan_gene_chr = gff_split[12]
        gmap_genome_name = gff_split[1]
        pan_genome_name = gff_split[13]
        gmap_stand_direction = gff_split[6]
        pan_gene_strand_direction = gff_split[18]
        #print(pan_genome_name)
        if gmap_chr == pan_gene_chr and gmap_genome_name == pan_genome_name and gmap_stand_direction == pan_gene_strand_direction:
            pan_gene_ID = gff_split[8].split(".")[0].split("=")[1]
            gmap_genome = gff_split[1]
            gmap_coordinate =  gmap_chr +":"+gff_split[3] + "-" + gff_split[3]
            gene_model_associate = gff_split[20].split(";")[0].split("=")[1]
            gene_info = (pan_gene_ID + "_" + gmap_genome + "," + gene_model_associate + ";" + gmap_coordinate)
            #print(gene_info)
            gff_intersect.append(gene_info)
            
            
with open('/Users/yinjieqiu/Desktop/pan_genome_nov 2/CDS_intersect_cutoff_90_all.txt', 'w') as f:
    for item in gff_intersect:
        f.write("%s\n" % item)
        
# NOTE, this file contains all gene models, which means gene ID from non-canonical transcripts are also included. To make it consistent, lines with only canonical transcripts are kept

# getting all canonical transcript ID: 
~/nam_pan_genome/NAM_annotation/canonical_transcript
cat *transcript.txt > all_canonical_transcript.txt

# a list of all gene model from CDS:
cut -d "," -f 2 CDS_intersect_cutoff_90_all.txt > model_list_to_filter.txt 
# convert p into to so can search with transcript ID
sed -i 's/_P0/_T0/g' model_list_to_filter.txt 
# 230548 lines


grep -Fwf ~/nam_pan_genome/NAM_annotation/canonical_transcript/all_canonical_transcript.txt model_list_to_filter.txt > canonical_model_gmap.txt

# 177249 lines are the canonical sequences  



# Final extraction: 
sed -i 's/0_P/0_T/g' CDS_intersect_cutoff_90_all.txt 

grep -Fwf canonical_model_gmap.txt CDS_intersect_cutoff_90_all.txt > Final_CDS_canonical_ID_replace_gmap9090_cutoff.txt

# the Final_CDS_canonical_ID_replace_gmap9090_cutoff.txt will be used to replace ID name before reshaping into R matrix 


# making dictionary for gff files to speed up the search 
# this will be the list made from above 
fill_in_info = {}
with open('/Users/yinjieqiu/Desktop/pan_genome_nov 2/Final_CDS_canonical_ID_replace_gmap9090_cutoff.txt') as all_info:
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
with open("/Users/yinjieqiu/Desktop/pan_genome_nov 2/gmap_info_for_fill_90.txt") as all_info:
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

with open('/Users/yinjieqiu/Desktop/pan_genome_nov 2/CDS_canonical.txt', 'w') as f:
    for item in final_list:
        f.write("%s\n" % item)


