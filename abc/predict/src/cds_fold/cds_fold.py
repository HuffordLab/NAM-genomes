from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description="Generate list of 0-fold and 4-fold sites from a fasta formatted reference genome and gff annotation file.")

parser.add_argument("reference", type=str, help="Fasta formatted reference file")

parser.add_argument("gff", type=str, help="version 3 gff file.")

parser.add_argument("-o", "--outfile", type=str, help="Name of the output file", required = True)

args = parser.parse_args()

###############
## FUNCTIONS ##
###############

def translate(seq):    
    seq = seq.upper()
    
    if len(seq) % 3 is not 0: 
        raise ValueError('input sequence not divisible by 3')
        
    if False in [s in ["A", "T", "G", "C", "N"] for s in seq]:
        raise ValueError(f'{seq} contains at least one non nucleotide or N character')
        
    table = { 
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W', 
    } 
    
    protein = "" 
    for i in range(0, len(seq), 3): 
        codon = seq[i:i + 3]
        try:
            protein += table[codon]
        except:
            protein += "-"
    return protein

   
#check all mutations at each position to determine fold
def condon(seq, pos):
    if len(seq) is not 3: 
        raise ValueError('input sequence should be 3 bps')
    if pos not in [1,2,3]:
        raise ValueError('pos argument must be 1, 2, or 3')
        
    if pos == 1: return [translate(b + seq[1:]) for b in ["A", "T", "G", "C"]] 
    if pos == 2: return [translate(seq[0] + b + seq[2]) for b in ["A", "T", "G", "C"]] 
    if pos == 3: return [translate(seq[0:2] + b) for b in ["A", "T", "G", "C"]]
    
#count potential number unique amino acids present at a locus
def n_fold(seq, gene_name):
    if len(seq) % 3 is not 0:
        print(f"skipping gene: {gene_name}. input sequence not divisible by 3")
        return None
    else:
        return [len(set(condon(seq[i:i+3], j))) for i in range(0, len(seq), 3) for j in [1,2,3]]

#quality check for issues with protein sequence
def build_seq(c_sequence, seqid, start, end, strand):
    if strand is "-":
        c_seq = fasta_dict[seqid].seq[start-1:end]
        return c_seq.reverse_complement() + c_sequence
    if strand is "+":
        c_seq = fasta_dict[seqid].seq[start-1:end]
        return c_sequence + c_seq

#quality check for issues with protein sequence
def protein_qc(seq, gene_name):

    if len(seq) % 3 is not 0:
        print(f"Warning! Skipping gene {gene_name}. Protein is not divisible by 3")
        return False

    if len(seq) < 3:
        print(f"Warning! Skipping gene {gene_name}. DNA sequence is too short.")
        return False

    protein = translate(seq)
    if len(protein) < 1:
        print(f"Warning! Skipping gene {gene_name}. Protein sequence is too short.")
        return False

    n_stops = protein.count("_")
    start_aa = protein[0]
    last_aa = protein[-1]

    if n_stops is not 1 or last_aa is not "_":
        print(f"Warning! Skipping gene: {gene_name}. Incorrect number and/or location of stop codons.")
        return False

    if start_aa is not "M":
        print(f"Warning! Skipping gene: {gene_name}. Does not start with M")
        return False

    if "-" in protein:
        print(f"Warning! Skipping gene: {gene_name}. Contains a non-amino acid three base pair condon")
        return False
    return True


#add to dictionary to track fold types at reference position
def append_fold(fold_dict, pos_key, fold):
    if pos_key not in fold_dict:
        fold_dict[pos_key] = []
        fold_dict[pos_key].append(fold)
    else:
        fold_dict[pos_key].append(fold)

#########################################################################
## READ REFERENCE GENOME SEQUENCE INTO A BIOPYTHON SEQUENCE DICITONARY ##
#########################################################################

fasta_dict = SeqIO.to_dict(SeqIO.parse(args.reference, "fasta"))


#################
## PROCESS GFF ##
#################

s = 0
smx = 1
old_name = ""
old_strand = ""
old_seqid = ""
c_sequence = ""
seq_pos = []
fold_dict = {}


with open(args.gff) as gff:
    for line in gff:
        if line[0] is not "#":
            seqid, source, s_type, start, end, score, strand, phase, attr = line.strip().split('\t')
            start = int(start)
            end = int(end)
            if s_type == "CDS":
                name = attr.split("ID=")[1].split(";")[0]
                if old_name == "": old_name = name
                if old_strand == "": old_strand = strand
                if old_seqid == "" : old_seqid = seqid
                if name == old_name:
                    c_sequence = build_seq(c_sequence, seqid, start, end, strand)
                    seq_pos += list(range(start, end+1))
                else:
                    gene_folds = n_fold(c_sequence, old_name)
                    prot_qc = protein_qc(c_sequence, old_name)
                    if prot_qc:
                        if gene_folds:
                            n_folds = len(gene_folds)
                            if old_strand == "-":
                                gene_folds = gene_folds[::-1]
                            for idx, fold in enumerate(gene_folds):
                                if fold == 1:
                                    append_fold(fold_dict, f"{old_seqid} {seq_pos[idx]}", 4)
                                if fold == 4:
                                    append_fold(fold_dict, f"{old_seqid} {seq_pos[idx]}", 0)
                    old_name = name
                    old_strand = strand
                    old_seqid = seqid
                    c_sequence = ""
                    c_sequence = build_seq(c_sequence, seqid, start, end, strand)
                    seq_pos = list(range(start, end+1))


outfile = open(args.outfile, 'w')
print(f"#chrom pos fold", file = outfile)
#print(f"chrom pos fold", file = zero_fold)
for k,v, in fold_dict.items():
    if len(set(v)) == 1:
        fold_one = list(set(v))[0]
        print(k, fold_one, file = outfile)

