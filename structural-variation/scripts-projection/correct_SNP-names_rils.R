#### arguments for command line ----

args <- commandArgs(trailingOnly = TRUE)

# # help
# if (all(length(args) == 1 & args == "-h" | args == "--help")) {
#   cat("
# Description: 
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
gbs.parents.file <- args[1]
gbs.rils.file <- args[2]

# gbs.parents.file <- "B73xB97/NAM_gbs-parents_SNPs.B73xB97.not-in-SVs.reseq-overlay.hmp.txt"
# gbs.rils.file <- "B73xB97/NAM_rils_SNPs.B73xB97.not-in-SVs.not-imputed.best-markers.hmp.txt"

#### libraries ----

if(!require("data.table")) install.packages("data.table")




#### load data ----

gbs.parents <- fread(gbs.parents.file, header = TRUE, data.table = FALSE)
gbs.rils <- fread(gbs.rils.file, header = TRUE, data.table = FALSE)



#### change SNP names ----

if (NROW(gbs.rils) == NROW(gbs.parents)) {
  if (all(gbs.rils[, "pos"] == gbs.parents[ , "pos"])) {
    cat("Changing marker names...\n")
    gbs.rils[ , 1] <- gbs.parents[, 1]
    cat("Done!\n")
  } else {
    cat("Changing marker names...\n")
    for (chr in unique(gbs.rils[, "chrom"])) {
      pos.rils <- gbs.rils[which(gbs.rils[, "chrom"] == chr), "pos"]
      gbs.rils[which(gbs.rils[, "chrom"] == chr), 1] <- gbs.parents[which(gbs.parents[, "chrom"] == chr & gbs.parents[, "pos"] %in% pos.rils), 1]
    }
    cat("Done!\n")
  }
} else {
  cat("Changing marker names...\n")
  for (chr in unique(gbs.rils[, "chrom"])) {
    pos.rils <- gbs.rils[which(gbs.rils[, "chrom"] == chr), "pos"]
    gbs.rils[which(gbs.rils[, "chrom"] == chr), 1] <- gbs.parents[which(gbs.parents[, "chrom"] == chr & gbs.parents[, "pos"] %in% pos.rils), 1]
  }
  cat("Done!\n")
}


# print new file
output.file <- gsub(pattern = ".hmp.txt", replacement = ".correct-marker-names.hmp.txt",
                    x = gbs.rils.file, fixed = TRUE)
cat("Writing ", output.file, "\n")
fwrite(gbs.rils, file = output.file, quote = FALSE, sep = "\t", na = NA, row.names = FALSE)
cat("Done!\n")
