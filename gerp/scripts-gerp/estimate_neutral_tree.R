library("rphast")
args <- commandArgs(trailingOnly = TRUE)
ref_root <- args[1]
msa_root <- args[2]
feat_file <- args[3]
tree <- args[4]
tree_output <- args[5]

chr_vector <- paste("chr", c(1:10), sep = "")

chr <- chr_vector[1]

ref_file <- paste(ref_root, chr, ".fa", sep = "")

msa <- paste(msa_root, chr, ".fa", sep="")

align4d <- read.msa(msa, refseq = ref_file, format = "FASTA", do.4d = TRUE, features= read.feat(feat_file))

for (chr in chr_vector[2:length(chr_vector)]){
  ref_file <- paste(ref_root, chr, ".fa", sep = "")
  msa <- paste(msa_root, chr, ".fa", sep="")
  align4d2 <- read.msa(msa, refseq = ref_file, format = "FASTA", do.4d = TRUE, features= read.feat(feat_file))
  align4d <- concat.msa(list(align4d, align4d2))
}

neutralMod <- phyloFit(align4d, tree=tree, subst.mod="REV")
sink(tree_output)
neutralMod$tree
sink()
