import os
import sys
import tempfile
import subprocess
import pandas as pd
import getopt
import re

#re added for gene name split - need to split on _ and . 

#Christine edited this version to be used in with the NAM genomes
#There are some issues with how the gene names/numbers are split
# Listing all of the functions
def usage():
    print("""\n
        This is the usage function:
        All_by_All_Blast.py
            -q or --query_gff          : Path to qeury gff file (genes only)
            -s or --subject_gff        : Path to subject gff file (genes only)
            -b or --blast_DB           : Path to subject blast DB of longest
                                         CDS transcript
            -l or --longest_CDS        : path to fasta file of subjects longest
                                         representative transcript sequence
            -o or --output_file        : output file name
            -n or --nucmer_coords      : Path to nucmer .coords file between
                                         subject and query
            -g or --subject_char_split : Text delimiter of split gene character

            -h or --help               : help command
        \n""")


try:
    opts, args = getopt.getopt(sys.argv[1:],
                               "q:s:b:l:o:n:g:h", ["query_gff=",
                                                   "subject_gff=",
                                                   "blast_DB=",
                                                   "longest_CDS=",
                                                   "output_file="
                                                   "nucmer_coords="
                                                   "subject_char_split=",
                                                   ]
                               )

except getopt.GetoptError as err:
    print(err)
    usage()
    sys.exit(2)

for opt, arg in opts:
    if opt in ("-q", "--query_gff"):
        query_gff = arg
    elif opt in ("-s", "--subject_gff"):
        Subject_gff_file = arg
    elif opt in ("-b", "--blast_DB"):
        DB = arg
    elif opt in ("-l", "--longest_CDS"):
        CDS = arg
    elif opt in ("-o", "--output_file"):
        output_file = arg
    elif opt in ("-n", "--nucmer_coords"):
        nucmer_coords = arg
    elif opt in ("-g", "--subject_char_split"):
        gene_char_split = arg
    elif opt in ("-h", "--help"):
        usage()
        sys.exit(2)
    else:
        assert False, "unhandled option"


def subprocess_cmd(command):
    subprocess.call(command,
                    stdout=subprocess.PIPE,
                    shell=True)


def cds_blast(gene, cds_fasta, DB):
    # make a temp file
    gene_file = tempfile.NamedTemporaryFile(suffix="_temp.txt",
                                            prefix="gene_",
                                            mode='r+',
                                            delete=False,
                                            )
    # get the basename
    base = os.path.basename(gene_file.name)
    base_split = os.path.splitext(base)[0]
    # we just needed the unique name so we can close it
    # gene_file.close()

    #  sed command to pull gene from gff
    subprocess_cmd("sed -n '/" +
                   gene +
                   "/,/^>/p' " +
                   cds_fasta +
                   " | head -n -1 >> " +
                   base_split +
                   'query.fa')

    # blast to arabidopsis gene db
    subprocess_cmd('blastn ' +
                   '-db ' +
                   DB +
                   ' ' +
                   '-query ' +
                   base_split +
                   'query.fa ' +
                   '-out ' +
                   gene_file.name +
                   ' ' +
                   '-task blastn ' +
                   '-evalue 1e-10 ' +
                   '-outfmt 6')

    # get the length of the CDS fasta
    cds_fasta_file = open(base_split+"query.fa", "r")
    cds_length = 0
    for line in cds_fasta_file:
        line = line.strip()
        if line.startswith(">"):
            continue
        else:
            cds_length += len(line)
    cds_fasta_file.close()

    # open the blast results and return a set of genes
    blast_results_file = open(gene_file.name, "r")
    blast_genes = []
    pid = []
    qstart = []
    qend = []
    e_val = []
    for line in blast_results_file:
        line = line.strip()
        fields = line.split()
        blast_genes.append(fields[1])
        pid.append(fields[2])
        qstart.append(fields[6])
        qend.append(fields[7])
        e_val.append(fields[10])
    # make it a unique list
    # unique_blast_genes = set(blast_genes)

    # Check if query file had anything
    if os.stat((base_split + 'query.fa')).st_size == 0:
        blast_genes = "CDS sequence not found"
        # remove old files
        subprocess_cmd('rm ' + base_split + 'query.fa')
        subprocess_cmd('rm ' + gene_file.name)
        gene_file.close()
        return(blast_genes)
    else:
        # remove old files
        subprocess_cmd('rm ' + base_split + 'query.fa')
        subprocess_cmd('rm ' + gene_file.name)
        gene_file.close()
        return (blast_genes, pid, qstart, qend, e_val, cds_length)


#This is my version of gene_grep where I'm just going to try to grep the gene name without the -w flag
def gene_grep(gene, gff):
    process = subprocess.check_output('grep ' +
                                      gene +
                                      ' ' +
                                      gff,
                                      shell=True
                                      ).decode('utf-8')
    process = process.strip()
    fields = process.split("\t")
    chrom = fields[0]
    start = fields[3]
    stop = fields[4]
    direction = fields[6]
    # pull out the gene name
    gene = fields[8].split(";")[0].split("=")[1]
    return(gene, chrom, start, stop, direction)


#Below is Meesh's version of gene_grep
#I'm removing the -w option because some of these transcripts are followed by _ which is considered a word constitutient character
#def gene_grep(gene, gff):
#    process = subprocess.check_output('grep -w ' +
#                                      gene +
#                                      ' ' +
#                                      gff,
#                                      shell=True
#                                      ).decode('utf-8')
#    process = process.strip()
#    fields = process.split("\t")
#    chrom = fields[0]
#    start = fields[3]
#    stop = fields[4]
#    direction = fields[6]
#    # pull out the gene name                                                                                                                                                                    
#    gene = fields[8].split(";")[0].split("=")[1]
#    return(gene, chrom, start, stop, direction)

# End of functions, start script
# Get gene names from Gff file
Query_gff_file = open(query_gff, "r")

# Open results file
OutFile = open(output_file, "w")
# open the coordinates dataframe
df = pd.read_csv(nucmer_coords,
                 sep="\t",
                 index_col=None,
                 low_memory=False
                 )

# chromosome naming can confuse python since it can start as int then go to str
df.TAGS = df.TAGS.astype(str)
df.TAGS2 = df.TAGS2.astype(str)

File_Header = print("Query_gene",
                    "Q_Chr",
                    "Q_start",
                    "Q_stop",
                    "Q_Syn_Chr",
                    "Q_Syn_start",
                    "Q_Syn_stop",
                    "S_Syn_Chr",
                    "S_Syn_start",
                    "S_Syn_stop",
                    "#_genes_in_syntentic",
                    "Sytentic_genes",
                    "#_adjacent_genes",
                    "adjacent_genes",
                    "gene_type",
                    "overlap_metadata",
                    "opposite_strands",
                    "tandem_dup_percentage",
                    "CD2S_low",
                    "CD2S_high",
                    sep="\t",
                    file=OutFile
                    )


for QueryGene in Query_gff_file:
    # make sure there is no hidden header
    if not QueryGene.lstrip().startswith('#'):
        Q_fields = QueryGene.split()
        Q_chr = Q_fields[0]
        Q_start = Q_fields[3]
        Q_stop = Q_fields[4]
        Q_size = int(Q_stop) - int(Q_start)
        Q_gene = Q_fields[8].split(";")[0].split("=")[1]
        # This is just for test cases
        # Q_chr = '4'
        # Q_start = 50276
        # Q_stop = 52356
        # Q_gene = 'Zm00001d027233'
        Q_gene_meta = [Q_gene, Q_chr, Q_start, Q_stop]

        # Call the function from above
        Blast = cds_blast(Q_gene, CDS, DB)

        # Set vars as NA unless syntenic region not found
        CD2S_low = "NA"
        CD2S_high = "NA"
        opposite_strands = "NA"
        t_dup_percent = "NA"

        # check to see if query gene was not found
        if "C" in Blast[0]:
            print(*Q_gene_meta,
                  "NA",
                  "NA",
                  "NA",
                  "NA",
                  "NA",
                  "NA",
                  "NA",
                  "NA",
                  "NA",
                  "NA",
                  "CDS_sequence_for_query_gene_not_found",
                  "NA",
                  "NA",
                  "NA",
                  "NA",
                  "NA",
                  sep="\t",
                  file=OutFile
                  )
            Blast_genes = []
            continue
        else:
            Blast_genes = Blast[0]
            Blast_pid = Blast[1]
            Blast_Qstart = Blast[2]
            Blast_Qstop = Blast[3]
            Blast_e_val = Blast[4]
            cds_length = int(Blast[5])
            cds_range = range(1, (cds_length + 1))
            # Get Chr, start and stop info for genes
            Blast_metadata = []
            # we don't want to loop through duplicate blasts so we use a set
            # we keep the lists to filter out multiple blasts later
            for idx, B_gene in enumerate(set(Blast_genes)):
                gene = B_gene.split("_")[0]
                B_gene_info = gene_grep(gene, Subject_gff_file)
                # multiple blast can be the same gene so find all  occurances
                idxs = [i for i, x in enumerate(Blast_genes) if x == B_gene]
                if len(idxs) == 1:
                    B_gene_info = (*B_gene_info,
                                   Blast_pid[Blast_genes.index(B_gene)],
                                   Blast_Qstart[Blast_genes.index(B_gene)],
                                   Blast_Qstop[Blast_genes.index(B_gene)],
                                   Blast_e_val[Blast_genes.index(B_gene)],
                                   )
                elif len(idxs) > 1:
                    largest_val_index = 0
                    largest_blast_size = 0
                    smallest_e_val = float(10000)
                    # Filter for multiple blasts to same gene by eval then len
                    for idx in idxs:
                        blast_size = (int(Blast_Qstop[idx]) -
                                      int(Blast_Qstart[idx]))
                        current_e_val = float(Blast_e_val[idx])
                        if current_e_val < smallest_e_val:
                            largest_val_index = idx
                            smallest_e_val = current_e_val
                        # in the off chance of equal e-vals
                        elif current_e_val == smallest_e_val:
                            if blast_size > largest_blast_size:
                                largest_blast_size = blast_size
                                largest_val_index = idx

                    B_gene_info = (*B_gene_info,
                                   Blast_pid[largest_val_index],
                                   Blast_Qstart[largest_val_index],
                                   Blast_Qstop[largest_val_index],
                                   Blast_e_val[largest_val_index],
                                   )
                Blast_metadata.append(list(B_gene_info))

            # find region of syntenty for the query gene
            try:
                Dframe_index = df.index[(df["TAGS"] == Q_chr) &
                                        (df['S1'] <= int(Q_stop))].tolist()[-1]
                Syn = df.iloc[(Dframe_index)].tolist()
                Two_sided = True

            # This pipeline will fail if the first gene in a Chr is not in a
                # Syntenic region, this is why we set this index exception
            # If this fails then we can only do a one sided one
            except IndexError as e:
                print("type error: ",
                      str(e),
                      "for gene",
                      Q_gene,
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
                        Syn_region_stop,
                        ]

            Syn_genes = []
            for result in Blast_metadata:
                gene = result[0]
                chrom = result[1]
                start = int(result[2])
                stop = int(result[3])
                if chrom == Syn_region_chrom:
                    if ((Syn_region_start)-500000) <= start <= \
                       ((Syn_region_stop) + 500000) and \
                       ((Syn_region_start)-500000) <= stop <= \
                       ((Syn_region_stop) + 500000):
                        Syn_genes.append(result)
            Syn_genes.sort(key=lambda x: x[0])

            if len(Syn_genes) > 1:
                adjacent_genes = []
                coordinates = []
                # print(Syn_genes)
                for idx, gene in enumerate(Syn_genes):
                    #coordinates.append(gene[0].split(gene_char_split)[1])
                    #above line is the issue (I think)
                    #print(gene[0])
                    gene_split = re.split('_|\.', gene[0])
                    #print(gene_split[0])
                    #above line allows me to split on gene name 'suffix' and then on divider character below
                    coordinates.append(gene_split[0].split(gene_char_split)[1])
                    #2020/07/15: with newest version of gene annotation the gene numbers are 10 apart, instead of 1
                    #So I changes the int(coordinates[idx]) - int(coordinates[idx-1]) to < 30 from < 3
                    if idx > 0:
                        if int(coordinates[idx]) - int(coordinates[idx-1]) < 110:
                            if Syn_genes[idx-1] not in adjacent_genes:
                                adjacent_genes.append(Syn_genes[idx-1])
                            adjacent_genes.append(Syn_genes[idx])

                if len(adjacent_genes) == 0:
                    gene_type = "non-adjacent_syntenic"
                    adjacent_genes_len = 0
                    overlap_metadata = "NA"
                else:
                    gene_type = "adjacent_genes_syntenic"
                    adjacent_genes_len = len(adjacent_genes)
                    ranges = []
                    stranded = []
                    s_cds_range = set(cds_range)
                    for gene in adjacent_genes:
                        blast_range = list(range(int(gene[6]),
                                           (1+int(gene[7]))))
                        ranges.append(blast_range)
                        stranded.append(gene[4])

                    # see if genes are on oppisite strands
                    if len(set(stranded)) == 2:
                        opposite_strands = True
                    else:
                        opposite_strands = False

                    if adjacent_genes_len == 2:
                        # use sets because we are useing buil in function
                        R1 = set(ranges[0])
                        R2 = set(ranges[1])
                        two_overlap = len(set.intersection(s_cds_range,
                                                           R1,
                                                           R2
                                                           ))
                        one_overlap = (len(set.intersection(s_cds_range, R1)) +
                                       len(set.intersection(s_cds_range, R2)) -
                                       (two_overlap*2
                                        ))
                        zero_overlap = cds_length - one_overlap - two_overlap
                        three_overlap = "NA"

                        # tandem dup percentage calc
                        t_dup_percent = (two_overlap) / (two_overlap +
                                                         one_overlap)
                        overlap_list = [str(cds_length),
                                        str(zero_overlap),
                                        str(one_overlap),
                                        str(two_overlap),
                                        str(three_overlap)]
                        overlap_metadata = ":".join(overlap_list)

                    elif adjacent_genes_len == 3:
                        R1 = set(ranges[0])
                        R2 = set(ranges[1])
                        R3 = set(ranges[2])
                        three_overlap = len(set.intersection(s_cds_range,
                                                             R1,
                                                             R2,
                                                             R3,))
                        two_overlap = (len(set.intersection(s_cds_range,
                                                            R1,
                                                            R2,
                                                            )) +
                                       len(set.intersection(s_cds_range,
                                                            R2,
                                                            R3)) +
                                       len(set.intersection(s_cds_range,
                                                            R1,
                                                            R3)) -
                                       (three_overlap*3))

                        one_overlap = (len(set.intersection(s_cds_range, R1)) +
                                       len(set.intersection(s_cds_range, R2)) +
                                       len(set.intersection(s_cds_range, R3)) -
                                       (three_overlap*3) - (two_overlap * 2))
                        zero_overlap = (cds_length -
                                        one_overlap -
                                        two_overlap -
                                        three_overlap)

                        # tandem dup percentage calc
                        t_dup_percent = (three_overlap) / (three_overlap +
                                                           two_overlap +
                                                           one_overlap)

                        overlap_list = [str(cds_length),
                                        str(zero_overlap),
                                        str(one_overlap),
                                        str(two_overlap),
                                        str(three_overlap)]
                        overlap_metadata = ":".join(overlap_list)

                    else:
                        overlap_list = [str(cds_length),
                                        "too many combinations",
                                        "NA",
                                        "NA"]
                        overlap_metadata = ":".join(overlap_list)

            elif len(Syn_genes) == 1:
                gene_type = "one_to_one_mapping"
                adjacent_genes = "NA"
                adjacent_genes_len = 0
                blast_range = range(int(Syn_genes[0][6]),
                                    (1+int(Syn_genes[0][7])))
                s_cds_range = set(cds_range)
                zero_overlap = (cds_length -
                                len(s_cds_range.intersection(blast_range)))
                one_overlap = len(s_cds_range.intersection(blast_range))
                two_overlap = "NA"
                three_overlap = "NA"
                overlap_list = [str(cds_length),
                                str(zero_overlap),
                                str(one_overlap),
                                str(two_overlap),
                                str(three_overlap)]
                overlap_metadata = ":".join(overlap_list)

            # Check to see if blast genes are near syntenic region
            else:
                # Set arbitrarily large values as a place holder
                CD2S_low = 1000000000000000
                CD2S_high = 1000000000000000
                # Check to see if you need to look on both sides of syn region
                if Two_sided:
                    for B_gene in Blast_metadata:
                        Before_candidate_list = []
                        After_candidate_list = []
                        if B_gene[1] == Syn_region_chrom:
                            before = Syn_region_start - int(B_gene[3])
                            after = int(B_gene[2]) - Syn_region_stop
                            if (before > 0) & (before < CD2S_low):
                                CD2S_low = before
                            elif (after > 0) & (after < CD2S_high):
                                CD2S_high = after
                            else:
                                continue
                else:
                    for index, B_gene in enumerate(Blast_metadata):
                        After_candidate_list = []
                        if B_gene[1] == Syn_region_chrom:
                            after = int(B_gene[2]) - Syn_region_stop
                            if (after > 0) & (after < CD2S_high):
                                CD2S_high = after
                            else:
                                continue
                if CD2S_low == 1000000000000000:
                    CD2S_low = "NA no CD2S low"
                if CD2S_high == 1000000000000000:
                    CD2S_high = "NA no CD2S high"

                gene_type = "non-matching_syntenic"
                adjacent_genes = "NA"
                adjacent_genes_len = 0
                overlap_metadata = "NA"

            # format the Syn gene lists
            P_Syn_genes = []
            for gene in Syn_genes:
                P_Syn_genes.append(",".join(map(str, gene)))

            P_adjacent_genes = []
            for gene in adjacent_genes:
                P_adjacent_genes.append(",".join(map(str, gene)))

            print(*Q_gene_meta,
                  *Syn_Meta,
                  len(Syn_genes),
                  ";".join(map(str, P_Syn_genes)),
                  adjacent_genes_len,
                  ";".join(map(str, P_adjacent_genes)),
                  gene_type,
                  overlap_metadata,
                  opposite_strands,
                  t_dup_percent,
                  CD2S_low,
                  CD2S_high,
                  sep="\t",
                  file=OutFile
                  )
            OutFile.flush()
# print(Q_gene, Q_start, Q_stop, Syn_genes)
OutFile.close()
