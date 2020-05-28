#!/usr/bin/python3

'''
created by Rafael Della Coletta
2019-09-30
'''


import argparse as ap


# initialize argument parser (pass user input from command line to script)
parser = ap.ArgumentParser(formatter_class=ap.RawDescriptionHelpFormatter,
                           description='''
description: this script reads a VCF file containing structural variant calls,
             and transform it into a hapmap file.

example: vcf2hapmap.py my_file.vcf my_results.hmp.txt''')
# add positional arguments
parser.add_argument("vcf_file", type=str,
                    help="VCF file with SV calls")
parser.add_argument("output_name", type=str,
                    help="name of the hapmap file")
# pass arguments into variables
args = parser.parse_args()
vcf_file = args.vcf_file
output_name = args.output_name


# open input file
infile = open(vcf_file, "r")

# open output file
outfile = open(output_name, "w")


# list with important info about SVs to extract from vcf file
vcf_fields = ["CHROM", "POS", "ID", "REF", "ALT",
              "QUAL", "FILTER", "INFO", "FORMAT"]
# create list to store indices of each info about snp/sv and inbreds in header
header_idx = []


# read file line by line
for line in infile:
    line = line.strip()
    # if line starts with ## -- skip
    if line[0:2] == "##":
        continue
    # if line starts with # -- save as header
    elif line[0:1] == "#":
        # save header
        header = line.split("\t")
        # remove # from #CHROM
        header[0] = header[0][1:]
        # get indices of columns with important information from header
        for field in vcf_fields:
            header_idx.append(header.index(field))
        # get indices for inbreds
        inbreds_list = header[9:]
        for inbred in inbreds_list:
            header_idx.append(header.index(inbred))
        # print(header_idx)
        # write output header
        print("rs", "alleles", "chrom", "pos", "strand", "assembly",
              "center", "protLSID", "assayLSID", "panel", "QCcode",
              "\t".join(inbreds_list), sep="\t", file=outfile)
    # otherwise, extract info from line
    else:
        line = line.split("\t")
        # get chrom number
        chr = line[header_idx[0]]
        # get length of sv
        sv_length = line[header_idx[7]].split(";SVLEN=")
        sv_length = sv_length[1].split(";")
        sv_length = abs(int(sv_length[0]))
        # get position in the middle of the SV
        sv_start = int(line[header_idx[1]])
        sv_end = sv_start + sv_length
        pos = round((sv_start + sv_end) / 2)
        # determine type of sv
        sv_type = line[header_idx[7]].split(";SVTYPE=")
        sv_type = sv_type[1].split(";")
        sv_type = sv_type[0].lower()
        # create id based on sv type and location
        id = sv_type + "." + chr + "." + str(sv_start) + "." + str(sv_end)
        # if sv is TRA, add TRA location in id
        if sv_type == "tra":
            # get location of where the TRA went
            tra_chr = line[header_idx[7]].split(";CHR2=")
            tra_chr = tra_chr[1].split(";")
            tra_chr = tra_chr[0]
            tra_pos = line[header_idx[7]].split(";END=")
            tra_pos = tra_pos[1].split(";")
            tra_pos = tra_pos[0]
            id = "tra." + tra_chr + "." + str(tra_pos)
        # parse each inbred line (based on its index on header)
        inbreds_geno = []
        for index in header_idx[9:]:
            # print(line[index])
            inbred_info = line[index].split(":")
            # if SV is 'T'here, assign genotype TT
            if inbred_info[0] == "1/1":
                genotype = "TT"
            # if missing info, assign genotype NN
            elif inbred_info[0] == "0/0" or inbred_info[0] == "0/1":
                genotype = "AA"
            # if SV is 'A'bsent or het, assign genotype AA
            else:
                genotype = "NN"
            inbreds_geno.append(genotype)
        # before writing output, format chrom according to hapmap format
        if chr[0:3] == "chr":
            chr = chr.split("chr")[1]
        if chr[0:4] == "scaf":
            chr = chr[0:4].upper() + chr[4:]
        # write output
        print(id, "A/T", chr, pos, "NA", "NA", "NA", "NA", "NA", "NA", "NA",
              "\t".join(inbreds_geno), sep="\t", file=outfile)


# close files
infile.close()
outfile.close()
