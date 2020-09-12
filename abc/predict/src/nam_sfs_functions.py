import pandas as pd
import re
import numpy as np
import gzip

#bed_df = pd.read_csv("data/ref/Zm-B73-REFERENCE-NAM-5.0_2Mb.bed", sep = "\t", names = ['chrome','start','end'])
#bed_df.query('@chr == chrome and (start <= @f_start <= end or start <= @f_end <= end)').index.tolist()
def openfile(filename):
    if filename.endswith(".gz"):
        return gzip.open(filename, "rt")
    else:
        return open(filename, "r")

#get dictionary sites to include. assume input is "chrome\tpos..."
def make_site_dict(site_file, fold_filter = None):
    site_dict = {}
    with open(site_file, "r") as sf:
        for line in sf:
            if fold_filter:
                chrom, pos, fold = line.strip().split()[0:3]
                #print(f"fold = {fold}")
            else:
                chrom, pos = line.strip().split()[0:2]
            chrpos = f"{chrom} {pos}"
            if fold_filter:
                if fold == fold_filter and chrpos not in site_dict:
                    site_dict[chrpos] = ""
            else:
                if chrpos not in site_dict:
                    site_dict[chrpos] = "" 
    return site_dict

#load windows
#load the bedfile made from bedtools makewindows into dictionary with keys chrom_start_end, 
    #values as empty sfs of specified size -- which is 26 for this project because we don't have invariants.
def build_sfsdf(bed_df, sfs_size):
    bed_idx = bed_df.index.tolist()
    sfs_dict = {}
    for idx in bed_idx:
        sfs_dict[idx] = [0] * sfs_size
    return sfs_dict


#NOT ADJUSTING FOR BED INDEXING WITH SV DATA
#which window of genome the variant falls into, INPUT IS BED INDEXING SO START IS VCF POS-1!
def query_bed_df(bed_df, q_chrome1, q_chrome2, q_start, q_end):
    q_start -= 1 #adjust for bed index
    if q_chrome1 == q_chrome2:
        idx = bed_df.query('@q_chrome1 == chrome and (start < @q_start <= end or start < @q_end <= end)').index.tolist()
    else: 
        q_end -= 1 #translocations aren't starts and ends, just single positions
        idx = bed_df.query('(@q_chrome1 == chrome and start < @q_start <= end) or (@q_chrome2 == chrome and  start < @q_end <= end)').index.tolist()
    return idx

#convert the position information to simple chr, pos pairs
#def sv_pos_parse(sv_pos):
#    chrom1, start, chrom2, end = re.split('\:|\-', sv_pos)
#    return chrom1, chrom2, start, end

#for new sv format, assumes input of list of whole line
def sv_pos_parse(sv_pos):
    rs, alleles, chrom2, end_or_mid = sv_pos[0:4]
    chrom2 = f"chr{chrom2}"
    #print(f"{rs}.{end_or_mid}")
    sv_type, chrom1, start, end = f"{rs}.{end_or_mid}".split(".")[0:4]
    return sv_type, chrom1, chrom2, start, end, alleles

#update folded_sfs
#input the current sfs and list of 0/1 formated genotypes
def sfs(sfs, geno_list, sfs_size, variant_fmt = "1/1"):

    #!!!! IF I SWITCH TO DOWNSAMPLING, DO IT HERE I THINK
    geno_sample = np.random.choice(geno_list, size = sfs_size, replace = False)

    #alt_category = sum([i.count("1/1") for i in geno_list]) 
    alt_category = sum([i.count(variant_fmt) for i in geno_sample]) 
    ref_category = len(sfs) - alt_category

    if alt_category > ref_category:
        sfs_category = ref_category
    else:
        sfs_category = alt_category

    #0th bin is singletons, python is 0-indexed
    if sfs_category > 0:
        sfs[sfs_category-1] += 1
    return sfs

#test
def sfs_test(sfs, geno_list):
    if geno_list >= 0:
        sfs[geno_list] += 1
    return sfs

#function to convert variant file line into list of "0/1" strings
#default is to allow NO heterzogous sites 
def line_parse(bed_df, sfs_dict, variant_line, 
               input_type, missing_prop, 
               heterozygous_prop, bi_allelic,
               i_count, site_dict = None, 
               sv_type = None):

    nucs = ["A", "T", "G", "C", "a", "t", "g", "c", "."]


    if input_type == "vcf":
  
        het_prop = float(variant_line.count("1/0") + variant_line.count("0/1")) / i_count
        miss_prop = float(variant_line.count("./.")) / i_count

        if het_prop <= heterozygous_prop and miss_prop <= missing_prop:

            line_list = variant_line.strip().split()
            variant_list = line_list[9:]
            if line_list[6] == "PASS":
                q_chrome =  line_list[0]
                q_end = int(line_list[1])
                q_start = int(q_end-1)
                ref, alt = line_list[2:4]
                #biallelic only
                
                site_pass = True 
                if bi_allelic:
                    site_pass = ref in nucs and alt in nucs
                #else:
                #    site_pass = True
                #!!! NOT WORKING!!!
                #print(site_pass)
                if site_pass:
                    if site_dict:
                        chr_pos = f"{q_chrome} {q_end}"
                        if chr_pos in site_dict:    
                            sfs_idx = query_bed_df(bed_df, q_chrome, q_chrome, q_start, q_end)
                            for idx in sfs_idx:
                                sfs_dict[idx] = sfs(sfs_dict[idx], variant_list, i_count)
                    else:
                        sfs_idx = query_bed_df(bed_df, q_chrome, q_chrome, q_start, q_end)
                        for idx in sfs_idx:
                            sfs_dict[idx] = sfs(sfs_dict[idx], variant_list, i_count)
            
    if input_type == "sv":

        het_prop = float(variant_line.count("AT") + variant_line.count("TA")) / i_count
        miss_prop = float(variant_line.count("NN")) / i_count

        #print(f"het and miss: {het_prop} {miss_prop}")

        if het_prop <= heterozygous_prop and miss_prop <= missing_prop:

            line_list = variant_line.strip().split()
            variant_list = line_list[11:]
            svtype, q_chrome1, q_chrome2, q_start, q_end, alleles = sv_pos_parse(line_list)
            sfs_idx = query_bed_df(bed_df, q_chrome1, q_chrome2, int(q_start), int(q_end))
            if len(alleles.split("/")) > 1:
                if sv_type:
                    if svtype in sv_type:
                        #print("LINE LIST IS:", line_list[0])
                        for idx in sfs_idx:
                            sfs_dict[idx] = sfs(sfs_dict[idx], variant_list, i_count, variant_fmt = "TT")
                else:
                        #print("LINE LIST IS:", line_list[0])
                        for idx in sfs_idx:
                            sfs_dict[idx] = sfs(sfs_dict[idx], variant_list, i_count, variant_fmt = "TT")

    return sfs_dict


def file_parse(in_file, input_type, bed_df, sfs_dict, i_count, bi_allelic = True,
               missing_prop = 0.0, heterozygous_prop = 0.0, site_dict = None, sv_type = None):
    with openfile(in_file) as f:
        for line in f:
            #print("line start is:", line[0:3])
            if "#" not in line and line[0:3] != "CHR":
                sfs_dict = line_parse(bed_df = bed_df, 
                           sfs_dict = sfs_dict,
                           variant_line = line,
                           input_type = input_type, 
                           bi_allelic = bi_allelic,
                           missing_prop = missing_prop,
                           heterozygous_prop = heterozygous_prop,
                           i_count = i_count, site_dict = site_dict, sv_type = sv_type)
    return sfs_dict

