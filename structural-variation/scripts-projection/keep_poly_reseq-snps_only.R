#### arguments for command line ----

args <- commandArgs(trailingOnly = TRUE)

# # help
# if (all(length(args) == 1 & args == "-h" | args == "--help")) {
#   cat("
# Description: this script merges SNPs and SVs hapmap files from usda parents to be used in
#              Tassel 5 when projecting SVs into RILs
#
# Usage: ")
#   quit()
# }

# make sure the correct number of arguments are used
# you should provide 2 arguments
if (length(args) != 2) {
  stop("incorrect number of arguments provided.

Usage:
       ")
}

# assign arguments to variables
cross <- args[1]
reseq.snps.file <- args[2]

# cross <- "B73xB97"
# reseq.snps.file <- "~/projects/sv_nams/data/tmp/B73xB97/NAM_parents-reseq_SNPs.B73xB97.not-in-SVs.hmp.txt"



#### libraries ----

if(!require("data.table")) install.packages("data.table")



#### keeping only polymorphic ----

# load data
reseq.snps <- fread(reseq.snps.file, header = TRUE, data.table = FALSE)

# get parents names
parent1 <- "B73"
parent2 <- toupper(unlist(strsplit(cross, "B73x"))[2])

# transform to all column names to upper case to avoid mismatch
colnames(reseq.snps)[12:NCOL(reseq.snps)]  <- toupper(colnames(reseq.snps)[12:NCOL(reseq.snps)])

# get parents column numbers in resequencing data
p1.col <- grep(parent1, colnames(reseq.snps))
p2.col <- grep(parent2, colnames(reseq.snps))

# get type of each marker
marker.type <- apply(X = reseq.snps[, c(p1.col, p2.col)],
                     MARGIN = 1, FUN = function(snp) {
                       
                       # get unique genotypes between parents
                       genotypes <- unique(snp)
                       genotypes <- genotypes[genotypes != "NN"]
                       
                       if (length(genotypes) == 0) {
                         
                         # if there is no genotype, snp is missing
                         return("missing")
                         
                       } else if (length(genotypes) == 1) {
                         
                         # if there is one genotype, it's monomorphic
                         # but distinguish if SNP is het
                         alleles <- unlist(strsplit(genotypes, split = ""))
                         if (alleles[1] == alleles[2]) {
                           return("mono")
                         } else {
                           return("het")
                         }
                         
                       } else {
                         
                         # if there are two genotypes, it's polymorphic
                         # but distiguish if one of the genotypes is het
                         p1.alleles <- unlist(strsplit(genotypes[1], split = ""))
                         p2.alleles <- unlist(strsplit(genotypes[2], split = ""))
                         if (p1.alleles[1] == p1.alleles[2] & p2.alleles[1] == p2.alleles[2]) {
                           return("poly")
                         } else {
                           return("het")
                         }
                         
                       }
                     })

# keep only homozygous polymorphic markers between parents
reseq.snps.poly <- reseq.snps[which(marker.type == "poly"), ]

# write results
reseq.snps.out <- gsub("not-in-SVs.hmp.txt", "poly.not-in-SVs.hmp.txt", reseq.snps.file)
fwrite(reseq.snps.poly, reseq.snps.out, quote = FALSE, sep = "\t",
       na = NA, row.names = FALSE)


