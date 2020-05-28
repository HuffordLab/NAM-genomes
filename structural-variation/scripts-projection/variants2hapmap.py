#!/usr/bin/python3

'''
created by Rafael Della Coletta
2020-02-04
'''


import argparse as ap


# initialize argument parser (pass user input from command line to script)
parser = ap.ArgumentParser(formatter_class=ap.RawDescriptionHelpFormatter,
                           description='''
description: this script reads a tab-delimeted file containing structural
             variant calls, and transform it into a hapmap file.

example: variants2hapmap.py my_file.txt my_results.sorted.hmp.txt''')
# add positional arguments
parser.add_argument("variant_file", type=str,
                    help="variant file with SV calls")
parser.add_argument("output_name", type=str,
                    help="name of the hapmap file")
# pass arguments into variables
args = parser.parse_args()
variant_file = args.variant_file
output_name = args.output_name


# open input file
infile = open(variant_file, "r")

# open output file
outfile = open(output_name, "w")

# get header information
header = infile.readline()
header = header.strip()
header = header.split("\t")
# get inbred list
inbreds_list = header[5:]

# write output header
print("rs", "alleles", "chrom", "pos", "strand", "assembly",
      "center", "protLSID", "assayLSID", "panel", "QCcode",
      "\t".join(inbreds_list), sep="\t", file=outfile)


# read file line by line
for line in infile:
    line = line.strip()
    line = line.split("\t")
    # get chrom number and start/end positions
    SV_location = line[0].split(":")
    chr = SV_location[0].split("-")[0]
    sv_start = int(SV_location[0].split("-")[1])
    sv_end = int(SV_location[1].split("-")[1])
    # get position in the middle of the SV
    pos = round((sv_start + sv_end) / 2)
    # # get length of sv
    # sv_length = abs(int(line[2]))
    # determine type of sv
    sv_type = line[4].lower()
    # create id based on sv type and location
    id = sv_type + "." + chr + "." + str(sv_start) + "." + str(sv_end)
    # if sv is TRA, add TRA location in id
    if sv_type == "tra":
        # get location of where the TRA went
        tra_chr = SV_location[1].split("-")[0]
        # correct id
        id = "tra." + tra_chr + "." + str(sv_end)
        # make sure position of TRA is the sv start, and not middle position
        pos = sv_start
    # parse each inbred line (based on its index on header)
    inbreds_geno = []
    for index in range(5, len(header)):
        # print(line[index])
        inbred_info = line[index]
        # if SV is 'T'here, assign genotype TT
        if inbred_info == "1/1":
            genotype = "TT"
        # if SV is 'A'bsent or het, assign genotype AA
        elif inbred_info == "0/0" or inbred_info == "0/1":
            genotype = "AA"
        # if missing info, assign genotype NN
        else:
            genotype = "NN"
        inbreds_geno.append(genotype)
    # before writing output, format chrom according to hapmap format
    if chr[0:3] == "chr":
        chr = chr.split("chr")[1]
    if chr[0:4] == "scaf":
        chr = chr[0:4].upper() + chr[4:]
    # check which alleles are present for that SV
    SV_alleles = "".join(inbreds_geno)
    if ("A" in SV_alleles) and ("T" in SV_alleles):
        alleles = "A/T"
    elif "A" in SV_alleles:
            alleles = "A"
    elif "T" in SV_alleles:
            alleles = "T"
    else:
            alleles = "N"
    # write output
    print(id, alleles, chr, pos, "NA", "NA", "NA", "NA", "NA", "NA", "NA",
          "\t".join(inbreds_geno), sep="\t", file=outfile)


# close files
infile.close()
outfile.close()
