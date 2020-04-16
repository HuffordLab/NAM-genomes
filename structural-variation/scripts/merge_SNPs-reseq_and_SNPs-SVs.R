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
# you should provide 5 arguments
if (length(args) != 5) {
  stop("incorrect number of arguments provided.

Usage:
       ")
}

# assign arguments to variables
snps.reseq <- args[1]
gbs.parents.file <- args[2]
gbs.rils.file <- args[3]
output.parents <- args[4]
output.rils <- args[5]


# snps.reseq <- "~/projects/sv_nams/data/tmp/B73xHp301/NAM_parents-reseq_SNPs.B73xHp301.not-in-SVs.hmp.txt"
# gbs.parents.file <- "~/projects/sv_nams/data/NAM_parents_SVs-SNPs.B73xHp301.sorted.hmp.txt"
# gbs.rils.file <- "~/projects/sv_nams/analysis/projection/NAM_rils_SVs-SNPs.B73xHp301.best-markers.projected.hmp.txt"
# output.parents <- "~/projects/sv_nams/data/NAM_parents_SNPs-reseq_and_SVs-SNPs.B73xHp301.hmp.txt"
# output.rils <- "~/projects/sv_nams/data/NAM_rils_SNPs-reseq_and_SVs-SNPs.B73xHp301.not-projected.hmp.txt"




#### libraries ----

if(!require("data.table")) install.packages("data.table")



#### function ----

MergeHapmaps <- function(snp_file, reseq_file, merge_RILs = FALSE, chr = NULL) {

  # load SNP and SV datasets
  snp.hmp <- fread(snp_file, header = TRUE, data.table = FALSE)
  sv.hmp <- fread(reseq_file, header = TRUE, data.table = FALSE)

  # correct B73 name in SV file
  colnames(sv.hmp)[grep("B73", colnames(sv.hmp))] <- "B73"
  
  # make sure columns are upper case to avoid mismatches
  colnames(snp.hmp)[12:NCOL(snp.hmp)] <- toupper(colnames(snp.hmp)[12:NCOL(snp.hmp)])
  colnames(sv.hmp)[12:NCOL(sv.hmp)] <- toupper(colnames(sv.hmp)[12:NCOL(sv.hmp)])
  
  # filter SVs that are missing in bith parents
  notB73.parent <- colnames(sv.hmp)[12:13][grep("B73", colnames(sv.hmp)[12:13], invert = TRUE)]
  sv.hmp <- sv.hmp[which(sv.hmp[, "B73"] != "NN" & sv.hmp[, notB73.parent] != "NN"), ]

  # if want to merge SNPs and SVs from parents...
  if (merge_RILs == FALSE) {

    # remove anything after "_" in SNP colnames
    for (i in 12:NCOL(snp.hmp)) {
      colnames(snp.hmp)[i] <- unlist(strsplit(colnames(snp.hmp)[i], split = "_"))[1]
    }

    # make sure columns have the same names and order in both files
    sv.hmp <- cbind(sv.hmp[, 1:11],
                    sv.hmp[, colnames(snp.hmp)[12:NCOL(snp.hmp)]])
    colnames(sv.hmp) <- colnames(snp.hmp)

    # merge hapmaps
    merged.hmp <- rbind(snp.hmp, sv.hmp)
  }

  # if want to create dataset for projection of SVs from parents to RILs...
  if (merge_RILs == TRUE) {
    
    # resequencing file is too big, so need to write files by chromosome
    if (chr != "scaffs") {
      snp.hmp <- subset(snp.hmp, chrom == chr)
      sv.hmp <- subset(sv.hmp, chrom == chr)
    } else {
      snp.hmp <- subset(snp.hmp, !chrom %in% 1:10)
      sv.hmp <- subset(sv.hmp, !chrom %in% 1:10)
    }
    
    # since i want to impute SVs from parents to RILs, I need to also include the SV marker names in
    # the RIL hapmap, but the actual genotypes of these markers are unknown (so far) and need to be
    # set to NA

    # create empty data frame with nrow = number of SVs, and ncol = number of RILs to be imputed
    rils.genos <- data.frame(matrix("NN", nrow = NROW(sv.hmp), ncol = length(12:NCOL(snp.hmp))))
    colnames(rils.genos) <- colnames(snp.hmp)[12:NCOL(snp.hmp)]

    # make sure values for alleles and genotypes are NAs before merging with RIL hmp
    rils.sv.hmp <- cbind(sv.hmp[, 1:11], rils.genos)
    rils.sv.hmp$alleles <- "NN"
    # also make sure the hapmap columns are the same
    colnames(rils.sv.hmp)[1:11] <- colnames(snp.hmp)[1:11]

    # merge hapmaps
    merged.hmp <- rbind(snp.hmp, rils.sv.hmp)
  }

  # sort by chromosome and position
  merged.hmp <- merged.hmp[order(merged.hmp$chrom, merged.hmp$pos),]

  return(merged.hmp)

}




#### merge hapmaps from parents ----

cat("Merging SNPs and SVs in parents...\n")
parents.merged <- MergeHapmaps(snp_file = gbs.parents.file, reseq_file = snps.reseq,
                               merge_RILs = FALSE)

fwrite(parents.merged, output.parents, quote = FALSE, sep = "\t", na = NA, row.names = FALSE)
cat("Done!\n")

#### merge hapmaps from RILs ----

for (chr in 1:10) {
  
  cat("Merging SNPs and SVs in RILs (", chr, ")...\n", sep = "")
  rils.merged <- MergeHapmaps(snp_file = gbs.rils.file, reseq_file = snps.reseq,
                              merge_RILs = TRUE, chr = chr)
  
  output.rils.chr <- gsub("not-projected.hmp.txt", paste0("chr-", chr, ".not-projected.hmp.txt"), output.rils)
  fwrite(rils.merged, output.rils.chr, quote = FALSE, sep = "\t", na = NA, row.names = FALSE)
  cat("Done!\n")
}

cat("Merging SNPs and SVs in RILs (scaffolds)...\n", sep = "")
rils.merged <- MergeHapmaps(snp_file = gbs.rils.file, reseq_file = snps.reseq,
                            merge_RILs = TRUE, chr = "scaffs")

output.rils.chr <- gsub("not-projected.hmp.txt", paste0("chr-scaffs.not-projected.hmp.txt"), output.rils)
fwrite(rils.merged, output.rils.chr, quote = FALSE, sep = "\t", na = NA, row.names = FALSE)
cat("Done!\n")



#### make sure parents and rils have the same markers ----

# if (NROW(parents.merged) > NROW(rils.merged)) {
#   cat("Making sure parents and rils have the same SNPs\n")
#   parents.merged.filtered <- parents.merged[which(parents.merged[, 1] %in% rils.merged[, 1]), ]
#   fwrite(parents.merged.filtered, output.parents, quote = FALSE, sep = "\t", na = NA, row.names = FALSE)
#   cat("Done!\n\n")
# }
