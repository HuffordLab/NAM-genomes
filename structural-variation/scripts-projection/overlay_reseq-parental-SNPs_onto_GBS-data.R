#### arguments for command line ----

args <- commandArgs(trailingOnly = TRUE)

# # help
# if (all(length(args) == 1 & args == "-h" | args == "--help")) {
#   cat("
# Description: this script adds parental SNPs from resequencing data into GBS, masking SNPs if
#              parental alleles from GBS and resequencing data disagree.
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
gbs.parents.file <- args[1]
reseq.file <- args[2]  # just need to specify chromosome 1 file (but other chromosomes should be in the same folder)
cross <- args[3]
output <- args[4]



#### libraries ----

if(!require("data.table")) install.packages("data.table")




#### load data ----

cat("Loading GBS data...\n")
gbs.data <- fread(gbs.parents.file, header = TRUE, data.table = FALSE)
cat("Done!\n\n")



#### overlay resequencing data onto gbs data ----

# make sure genotypes are in uppercase
colnames(gbs.data)[12:NCOL(gbs.data)] <- toupper(colnames(gbs.data)[12:NCOL(gbs.data)])

# get parent names
p1.gbs.col <- grep("B73_", colnames(gbs.data))
p1.gbs <- colnames(gbs.data)[p1.gbs.col]

p2.gbs <- toupper(unlist(strsplit(cross, split = "B73x"))[2])
p2.gbs.col <- grep(paste0(p2.gbs, "_"), colnames(gbs.data))
p2.gbs <- colnames(gbs.data)[p2.gbs.col]

# filter gbs data to have only the parents of cross being analyzed
gbs.parents <- gbs.data[, c(1:11, p1.gbs.col, p2.gbs.col)]

# create empty df to store results
gbs.parents.overlay <- data.frame(matrix(nrow = 0, ncol = NCOL(gbs.parents)))
colnames(gbs.parents.overlay) <- colnames(gbs.parents)

for (chr in unique(gbs.data[, "chrom"])) {

  cat("Overlaying SNPs from resequencing data on chr", chr, "...\n", sep = "")
  # subset gbs data by chromosome
  gbs.parents.chr <- gbs.parents[which(gbs.parents[, "chrom"] == chr), ]

  # load reseq file chromosome by chromosome
  if (chr %in% as.character(1:10)) {
    reseq.file.chr <- gsub(pattern = "chr[0-9]+", replacement = paste0("chr", chr), x = reseq.file,
                           perl = TRUE)
  } else {
    reseq.file.chr <- gsub(pattern = "chr[0-9]+", replacement = "scaffs", x = reseq.file,
                           perl = TRUE)
  }
  reseq.data.chr <- fread(reseq.file.chr, header = TRUE, data.table = FALSE)
  reseq.data.chr <- reseq.data.chr[which(reseq.data.chr[, "chrom"] == chr), ]
  reseq.data.chr <- reseq.data.chr[which(reseq.data.chr[, "pos"] %in% gbs.parents.chr[, "pos"]), ]
  # also make sure that gbs data has the snps as reseq now
  gbs.parents.chr <- gbs.parents.chr[which(gbs.parents.chr[, "pos"] %in% reseq.data.chr[, "pos"]), ]

  # filter reseq data to have only the parents of cross being analyzed
  # make sure they are all uppercase
  p1.reseq <- "B73"
  p2.reseq <- toupper(unlist(strsplit(cross, split = "B73x"))[2])
  colnames(reseq.data.chr)[12:NCOL(reseq.data.chr)] <- toupper(colnames(reseq.data.chr)[12:NCOL(reseq.data.chr)])
  reseq.data.chr <- cbind(reseq.data.chr[, 1:11], reseq.data.chr[, c(p1.reseq, p2.reseq)])

  # transform gbs SNP names into names of resequencing
  if (all(reseq.data.chr[, "pos"] == gbs.parents.chr[ , "pos"])) {
    gbs.parents.chr[ , 1] <- reseq.data.chr[, 1]
  }

  # get indices of SNP calls that disagree between reseq and gbs
  p1.disagree <- which(gbs.parents.chr[, p1.gbs] != reseq.data.chr[, p1.reseq])
  p2.disagree <- which(gbs.parents.chr[, p2.gbs] != reseq.data.chr[, p2.reseq])

  cat("  Parent ", p1.gbs, "...", sep = "")

  # overlay reseq data on missing gbs data on parent 1
  for (i in p1.disagree) {
    if (gbs.parents.chr[i, p1.gbs] == "NN") {
      # if missing on gbs data, use allele from reseq
      gbs.parents.chr[i, p1.gbs] <- reseq.data.chr[i, p1.reseq]
    } else if (gbs.parents.chr[i, p1.gbs] != "NN" & reseq.data.chr[i, p1.reseq] == "NN") {
      # if not missing on gbs data AND missing on reseq, keep allele from gbs
      gbs.parents.chr[i, p1.gbs] <- gbs.parents.chr[i, p1.gbs]
    } else {
      # if not missing on gbs data and not missing on reseq, transform disagreement as "NN"
      gbs.parents.chr[i, p1.gbs] <- "NN"
    }
  }

  cat(" done!\n  Parent ", p2.gbs, "...", sep = "")

  # repeat for parent 2
  for (i in p2.disagree) {
    if (gbs.parents.chr[i, p2.gbs] == "NN") {
      # if missing on gbs data, use allele from reseq
      gbs.parents.chr[i, p2.gbs] <- reseq.data.chr[i, p2.reseq]
    } else if (gbs.parents.chr[i, p2.gbs] != "NN" & reseq.data.chr[i, p2.reseq] == "NN") {
      # if not missing on gbs data AND missing on reseq, keep allele from gbs
      gbs.parents.chr[i, p2.gbs] <- gbs.parents.chr[i, p2.gbs]
    } else {
      # if not missing on gbs data and not missing on reseq, transform disagreement as "NN"
      gbs.parents.chr[i, p2.gbs] <- "NN"
    }
  }

  # cross B73xTzi8 doesn't have gbs data on parent Tzi8. Thus, I'll just use the resequencing data
  if (cross == "B73xTzi8") {
    gbs.parents.chr <- cbind(gbs.parents.chr, TZI8 = reseq.data.chr[, p2.reseq])
  }


  cat(" done!\n")

  gbs.parents.overlay <- rbind(gbs.parents.overlay, gbs.parents.chr)
}

# write output
fwrite(gbs.parents.overlay, output, quote = FALSE, sep = "\t", na = NA, row.names = FALSE)
