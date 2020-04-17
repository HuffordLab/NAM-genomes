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
chr <- args[1]
folder.after.proj <- args[2]



# chr <- 10
# folder.after.proj <- "~/projects/sv_nams/analysis/reseq_snps_projection"




#### libraries ----

if(!require("data.table")) install.packages("data.table")


# get list with all NAM families
cross.list <- system("ls -d ~/projects/sv_nams/data/GBS-output/tmp/B73* | xargs -n 1 basename",
                     intern = TRUE)

# get a list with all sv positions per chromosome for all crosses
all.proj.sv.pos <- list()

cat("Getting list with SV positions for all populations...\n")

for (cross in cross.list) {
  
  cat("  ", cross, "\n")
  
  # get parents names
  parent1 <- "B73"
  parent2 <- toupper(unlist(strsplit(cross, "B73x"))[2])
  
  # load hapmap after projection
  filename.after.proj <- list.files(path = folder.after.proj,
                                    pattern = paste0(cross, ".only-reseq-snps.chr-", chr,".projected.hmp.txt"),
                                    full.names = TRUE)
  hmp.after <- fread(filename.after.proj, header = TRUE, data.table = FALSE)
  # filter hmp by svs and select only chromosome and positions
  hmp.after <- hmp.after[, c(1, 3, 4)]
  
  if (!as.character(chr) %in% names(all.proj.sv.pos)) {
    all.proj.sv.pos[[as.character(chr)]] <- hmp.after[, "pos"]
  } else {
    all.proj.sv.pos[[as.character(chr)]] <- append(all.proj.sv.pos[[as.character(chr)]], hmp.after[, "pos"])
    all.proj.sv.pos[[as.character(chr)]] <- sort(unique(all.proj.sv.pos[[as.character(chr)]]))
  }
  
}

cat("Done!\n\n")

# create an empty hmp file with hmp columns
all.svs.hmp <- data.frame(matrix("NA", nrow = length(all.proj.sv.pos[[as.character(chr)]]),
                                 ncol = 11), stringsAsFactors = FALSE)
colnames(all.svs.hmp) <- c("rs#", "alleles", "chrom", "pos", "strand", "assembly#", "center", 
                           "protLSID", "assayLSID", "panelLSID", "QCcode" )
# fix first columns
all.svs.hmp[, "rs#"] <- paste0("S", chr, "_", all.proj.sv.pos[[as.character(chr)]])
all.svs.hmp[, "alleles"] <- "N"
all.svs.hmp[, "chrom"] <- chr
all.svs.hmp[, "pos"] <- all.proj.sv.pos[[as.character(chr)]]
all.svs.hmp[, "strand"] <- "+"

# write file
out.filename <- paste0(folder.after.proj, "/NAM_rils_SNPs-reseq_and_SVs-SNPs.all-reseq-snps.empty.chr-", chr, ".hmp.txt")
fwrite(all.svs.hmp, file = out.filename, quote = FALSE, sep = "\t", na = "NA", row.names = FALSE)


