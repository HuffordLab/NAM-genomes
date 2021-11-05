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
# you should provide 4 arguments
if (length(args) != 4) {
  stop("incorrect number of arguments provided.

Usage:
       ")
}

# assign arguments to variables
cross <- args[1]
chr <- args[2]
folder.not.proj <- args[3]
folder.after.proj <- args[4]


# cross <- "B73xB97"
# chr <- 10
# folder.not.proj <- "~/projects/sv_nams/data/tmp"
# folder.after.proj <- "~/projects/sv_nams/analysis/reseq_snps_projection2"




#### libraries ----

if(!require("data.table")) install.packages("data.table")




#### add monomorphic snps by chromosome ----

cat(cross, "- chr", chr, "\n")

# get parents names
parent1 <- "B73"
parent2 <- toupper(unlist(strsplit(cross, "B73x"))[2])

# get sv filename for that cross# get list with all SV filenames
sv.filename <- list.files(folder.not.proj,
                          pattern = paste0("NAM_parents-reseq_SNPs.", cross, ".not-in-SVs.hmp.txt"),
                          recursive = TRUE,
                          full.names = TRUE)
sv.filename <- sv.filename[grep(".gz", sv.filename, invert = TRUE)]
# open file with SV positions
sv.hmp.cross <- fread(sv.filename, header = TRUE, data.table = FALSE)
# filter by chromosome
sv.hmp.cross <- subset(sv.hmp.cross, chrom == chr)
# change parent columns to upper case
colnames(sv.hmp.cross)[12:NCOL(sv.hmp.cross)] <- toupper(colnames(sv.hmp.cross)[12:NCOL(sv.hmp.cross)])
# rename B73 parent
colnames(sv.hmp.cross)[grep("B73", colnames(sv.hmp.cross))] <- "B73"

# get type of each marker
marker.type <- apply(X = sv.hmp.cross[, c(parent1, parent2)],
                     MARGIN = 1, FUN = function(snp) {

                       # get unique genotypes between parents
                       genotypes <- unique(snp)

                       if (any(grepl("NN", genotypes))) {

                         # if there is NN, consider SNP is missing
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

# keep only homozygous monomorphic markers between parents
sv.hmp.cross.mono <- sv.hmp.cross[which(marker.type == "mono"), ]
# sv.hmp.cross.missing.het <- sv.hmp.cross[which(marker.type == "missing" | marker.type == "het"), ]

# after projection
filename.after.proj <- list.files(path = folder.after.proj,
                                  pattern = paste0("NAM_rils_SNPs-only.", cross),
                                  full.names = TRUE)
filename.after.proj <- filename.after.proj[grep("poly.projected.hmp.txt", filename.after.proj)]
hmp.after <- fread(filename.after.proj, header = TRUE, data.table = FALSE)
# filter by chromosome
hmp.after <- subset(hmp.after, chrom == chr)
# keep only reseq SNPs
hmp.after <- hmp.after[grep(paste0("^S", chr, "_"), hmp.after[, 1], perl = TRUE), ]
# remove duplicates
hmp.after <- hmp.after[!duplicated(hmp.after[, 1]), ]


# replicate columns of monomorphic parental SNPs to the have the number of RILs after projection
reps <- length(12:NCOL(hmp.after)) - 2
sv.hmp.cross.mono <- cbind(sv.hmp.cross.mono, replicate(reps, sv.hmp.cross.mono$B73))
# make sure columns have the same name as rils
colnames(sv.hmp.cross.mono) <- colnames(hmp.after)

# merge mono and polymorphic SNPs
merged.hmp <- rbind(sv.hmp.cross.mono, hmp.after)
# sort by chromosome and position
merged.hmp <- merged.hmp[order(merged.hmp$pos), ]

# keep only RILs
ril.columns <- grep("Z[0-9][0-9][0-9]E", colnames(merged.hmp), perl = TRUE)
merged.hmp <- cbind(merged.hmp[, 1:11], merged.hmp[, ril.columns])

# write merged files
outfile <- gsub("poly.projected.hmp.txt", paste0("only-reseq-snps.chr-", chr,".projected.hmp.txt"), filename.after.proj)
fwrite(merged.hmp, outfile, quote = FALSE, sep = "\t", na = NA, row.names = FALSE)
