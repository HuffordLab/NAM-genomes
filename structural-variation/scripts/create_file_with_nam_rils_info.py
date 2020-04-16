#!/usr/bin/python3

'''
version 1.0

created by Rafael Della Coletta
2019-10-03
'''


import argparse as ap


# initialize argument parser (pass user input from command line to script)
parser = ap.ArgumentParser(formatter_class=ap.RawDescriptionHelpFormatter,
                           description='''
description: this script creates a tab-delimited file with information about
             NAM crosses and their respective RIL populations.''')
# add positional arguments
parser.add_argument("file_with_info_about_NAMs", type=str,
                    help="hapmap file with SNPs calls")
# pass arguments into variables
args = parser.parse_args()
info_file = args.file_with_info_about_NAMs

# cross information based on http://maizecoop.cropsci.uiuc.edu/nam-rils.php
nam_dict = {"B73xB97":"Z001", "B73xCML103":"Z002", "B73xCML228":"Z003",
            "B73xCML247":"Z004", "B73xCML277":"Z005", "B73xCML322":"Z006",
            "B73xCML333":"Z007", "B73xCML52":"Z008", "B73xCML69":"Z009",
            "B73xHp301":"Z010",  "B73xIl14H":"Z011", "B73xKi11":"Z012",
            "B73xKi3":"Z013", "B73xKy21":"Z014", "B73xM162W":"Z015",
            "B73xM37W":"Z016", "B73xMo18W":"Z018", "B73xMS71":"Z019",
            "B73xNC350":"Z020", "B73xNC358":"Z021", "B73xOh43":"Z022",
            "B73xOh7B":"Z023", "B73xP39":"Z024", "B73xTx303":"Z025",
            "B73xTzi8":"Z026"}

# read file with list of NAM parents and RILs
with open(info_file, "r") as nam_file:
    # read line containing info about nam lines
    nam_lines = nam_file.readline()

# split lines by comma
nam_lines = nam_lines.strip()
nam_lines = nam_lines.split(",")

# get list of parents
parents_list = [nam for nam in nam_lines if nam.find("_") > -1]

# start a dictionary to save RILs for each cross
rils_from_cross = {}

# open output file (it's going to overwrite the input file)
outfile = open(info_file, "w")
# print header
print("cross", "parents", "rils", sep="\t", file=outfile)

for cross in nam_dict.keys():
    # get rils based on IDs of each cross
    id = nam_dict[cross].split(",")[0]
    rils_from_cross[cross] = [nam for nam in nam_lines if nam[0:4] == id]
    # get B73 ID as in the GBS data
    B73_parent = cross.split("x", maxsplit=1)[0]
    B73_id = [id for id in parents_list if id.find(B73_parent + "_") > -1]
    B73_id = "".join(B73_id)
    # get other NAM parent ID as in the GBS data
    nam_parent = cross.split("x", maxsplit=1)[1]
    nam_id = [id for id in parents_list if id.find(nam_parent + "_") > -1]
    if len(nam_id) == 0:
        # looks like Tzi8 inbred was not genotyped with GBS, so in this case i
        # will just add its name to complete the table
        nam_id = nam_parent
    else:
        nam_id = "".join(nam_id)
    # print output
    print(cross, B73_id,  nam_id, ",".join(rils_from_cross[cross]),
          sep="\t", file=outfile)

# close file
outfile.close()
