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
sv.file <- args[1]
gbs.parents.file <- args[2]
gbs.rils.file <- args[3]
output.parents <- args[4]
output.rils <- args[5]


# sv.file <- "~/projects/sv_nams/data/NAM_founders_SVs_B73xTzi8.hmp.txt"
# gbs.parents.file <- "B73xTzi8/NAM_gbs-parents_SNPs.B73xTzi8.not-in-SVs.reseq-overlay.hmp.txt"
# gbs.rils.file <- "B73xTzi8/NAM_rils_SNPs.B73xTzi8.not-in-SVs.not-imputed.best-markers.correct-marker-names.hmp.txt"
# output.parents <- "~/projects/sv_nams/data/NAM_parents_SVs-SNPs.B73xTzi8.hmp.txt"
# output.rils <- "~/projects/sv_nams/data/NAM_rils_SVs-SNPs.B73xTzi8.best-markers.not-projected.hmp.txt"


      
#### libraries ----

if(!require("data.table")) install.packages("data.table")



#### function ----

MergeHapmaps <- function(snp_file, parental_sv_file, merge_RILs = FALSE) {
  
  # load SNP and SV datasets
  snp.hmp <- fread(snp_file, header = TRUE, data.table = FALSE)
  sv.hmp <- fread(parental_sv_file, header = TRUE, data.table = FALSE)
  
  # remove ".bam" from SV colnames and make it all upper case letters
  colnames(sv.hmp)[12:NCOL(sv.hmp)] <- toupper(gsub(".bam", "", colnames(sv.hmp)[12:NCOL(sv.hmp)], fixed = TRUE))
  # correct B73 name in SV file
  colnames(sv.hmp)[grep("B73", colnames(sv.hmp))] <- "B73"
  
  # filter SVs that are missing in non-B73 parent
  notB73.parent <- colnames(sv.hmp)[12:13][grep("B73", colnames(sv.hmp)[12:13], invert = TRUE)]
  sv.hmp <- sv.hmp[which(sv.hmp[, notB73.parent] != "NN"), ]
  
  # B73 should not have any SV because all SVs were called against B73 ref genome
  # B73Ab10 was used as a negative control, but there is a lot of missing data
  # thus I will convert "NN" to "AA", and leave all non-missing calls as they were called
  sv.hmp[which(sv.hmp[, "B73"] == "NN"), "B73"]  <- "AA"

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
parents.merged <- MergeHapmaps(snp_file = gbs.parents.file, parental_sv_file = sv.file,
                               merge_RILs = FALSE)

fwrite(parents.merged, output.parents, quote = FALSE, sep = "\t", na = NA, row.names = FALSE)
cat("Done!\n")

#### merge hapmaps from RILs ----

cat("Merging SNPs and SVs in RILs\n")
rils.merged <- MergeHapmaps(snp_file = gbs.rils.file, parental_sv_file = sv.file,
                            merge_RILs = TRUE)

fwrite(rils.merged, output.rils, quote = FALSE, sep = "\t", na = NA, row.names = FALSE)
cat("Done!\n")



#### make sure parents and rils have the same markers ----

if (NROW(parents.merged) > NROW(rils.merged)) {
  cat("Making sure parents and rils have the same SNPs\n")
  parents.merged.filtered <- parents.merged[which(parents.merged[, 1] %in% rils.merged[, 1]), ]
  fwrite(parents.merged.filtered, output.parents, quote = FALSE, sep = "\t", na = NA, row.names = FALSE)
  cat("Done!\n\n")
}
