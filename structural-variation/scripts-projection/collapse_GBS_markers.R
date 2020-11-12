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
  
Usage: Rscript collapse_GBS_data [nam_file] [cross]
       ")
}

# assign arguments to variables
nam.file <- args[1]
cross <- args[2]

# nam.file <- "data/NAM_rils_SNPs.B73xB97.chr1.not-in-SVs.hmp.txt"
# cross <- "B73xB97"


#### libraries ----

if(!require("data.table")) install.packages("data.table")
if(!require("doParallel")) install.packages("doParallel")




#### load data ----

cat("Loading ", cross, "...\n", sep = "")
nam.hmp <- fread(nam.file, header = TRUE, data.table = FALSE)



#### check missing data ----

# get only few lines to test
p2 <- unlist(strsplit(cross, split = "B73x"))[2]

# get column number for parents, hybrids and RILs
col.number.b73 <- grep("B73_", colnames(nam.hmp))
col.number.p2 <- grep(paste0(p2, "_"), colnames(nam.hmp))
col.number.hyb <- grep(cross, colnames(nam.hmp))
col.number.rils <- grep("Z0[0-9][0-9]", colnames(nam.hmp))

# remove SNPs with b73 missing data
hmp.not.miss <- nam.hmp[which(nam.hmp[, col.number.b73] != "NN"), ]

cat(round((sum(nam.hmp[, col.number.b73] == "NN") / NROW(nam.hmp) * 100), digits = 2),
    "% of SNPs are missing on B73\n", sep="")
cat(round((sum(nam.hmp[, col.number.p2] == "NN") / NROW(nam.hmp) * 100), digits = 2),
    "% of SNPs are missing on parent 2\n", sep="")
cat(round((sum(nam.hmp[, col.number.hyb] == "NN") / NROW(nam.hmp) * 100), digits = 2),
    "% of SNPs are missing on hybrid\n", sep="")

NN.count.total <- c()
for (i in 16:NCOL(nam.hmp)) {
  NN.count <- sum(nam.hmp[, i] == "NN")
  NN.count.total <- append(NN.count.total, NN.count)
  # cat(colnames(nam.hmp)[i], "has", NN.count, "missing data\n")
}
cat("On average, ", round((mean(NN.count.total) / NROW(nam.hmp) * 100), digits = 2),
    "% of SNPs are missing on RILs\n\n", sep="")



#### collapse duplicated SNPs ----

nam.hmp.collapsed <- data.frame(matrix(nrow = 0, ncol = NCOL(nam.hmp)))
colnames(nam.hmp.collapsed) <- colnames(nam.hmp)


# find and collapse duplicates per chromosome
for (chr in unique(nam.hmp[, "chrom"])) {
  
  cat("Collapsing duplicated SNPs in chromosome ", chr, "...\n", sep = "")
  nam.hmp.chr <- nam.hmp[which(nam.hmp[, "chrom"] == chr), ]
  duplicates <- nam.hmp.chr[, "pos"][duplicated(nam.hmp.chr[, "pos"])]
  # since some rows can be repeated more than 2 times, need to get unique duplicated values
  duplicates <- unique(duplicates)
  
  # get number of cores to run in parallel
  num.cores <- detectCores()
  # get collapsed rows
  collapsed.rows <- mclapply(duplicates, FUN = function(dup.pos, chr.df) {
    
    # get all rows that have that position
    chr.df.dup <- chr.df[which(chr.df[, "pos"] == dup.pos), ]
    
    # collapse genotype calls for all lines
    collapsed.SNPs <- sapply(X = chr.df.dup[, 12:NCOL(chr.df.dup)],
                             FUN = function(geno) {
                               # if there's only one unique genotype called in all duplicates
                               if (length(unique(geno)) == 1) {
                                 # return that genotype
                                 return(unique(geno))
                                 
                                 # if there's more than one genotype call
                               } else {
                                 # check how many other non-missing genotypes there are
                                 not.NN.geno <- geno[which(geno != "NN")]
                                 # if there's only one unique genotype other than NN
                                 if (length(unique(not.NN.geno)) == 1) {
                                   # return that genotype
                                   return(unique(not.NN.geno))
                                 } else {
                                   # if there's more than one, it means there's a disagreement between calls,
                                   # so return "NN"
                                   return("NN")
                                 }
                               }
                             })
    
    # save information about SNPs
    SNP.info <- chr.df.dup[1, 1:11]
    
    filtered.row <- cbind(SNP.info, as.list(collapsed.SNPs), stringsAsFactors = FALSE)
    
    # return row of dataframe
    return(filtered.row)
    
  }, nam.hmp.chr, mc.cores = num.cores)
  
  # transform lists into dataframe
  collapsed.df <- do.call(rbind, collapsed.rows)
  # remove all duplicated rows from main df
  nam.hmp.chr.no.dups <- nam.hmp.chr[which(!nam.hmp.chr[, "pos"] %in% duplicates), ]
  # merge collapsed.df to main df without duplicates
  nam.hmp.chr.filtered <- rbind(nam.hmp.chr.no.dups, collapsed.df)
  # sort dataframeby chr and position
  nam.hmp.chr.filtered <- nam.hmp.chr.filtered[order(nam.hmp.chr.filtered[, "chrom"],
                                                     nam.hmp.chr.filtered[, "pos"]), ]
  
  # append to df with other chromosomes
  nam.hmp.collapsed <- rbind(nam.hmp.collapsed, nam.hmp.chr.filtered)
  
  cat("Done!\n\n")
}

# output name
outfile <- unlist(strsplit(nam.file, split = ".hmp.txt", fixed = TRUE))
# write table
fwrite(x = nam.hmp.collapsed, file = paste0(outfile, ".collapsed.hmp.txt"), na = NA,
       quote = FALSE, sep = "\t")

