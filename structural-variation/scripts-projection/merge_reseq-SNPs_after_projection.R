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



cat("Creating file with SVs only for each cross...\n")

for (cross in cross.list) {
  
  cat("  ", cross, "\n")
  
  # load hapmap after projection
  filename.after.proj <- list.files(path = folder.after.proj,
                                    pattern = paste0(cross, ".only-reseq-snps.chr-", chr,".projected.hmp.txt"),
                                    full.names = TRUE)
  hmp.after <- fread(filename.after.proj, header = TRUE, data.table = FALSE)
  
  # read empty hmp file with snps from all crosses
  all.svs.empty.hmp.file <- list.files(path = folder.after.proj,
                                       pattern = paste0("all-reseq-snps.empty.chr-", chr,".hmp.txt"),
                                       full.names = TRUE)
  all.svs.hmp <- fread(all.svs.empty.hmp.file, header = TRUE, data.table = FALSE)
  
  # View(all.svs.hmp[1:10, ])
  
  # cbind NN matrix with column number = number of genotypes in a cross
  all.svs.hmp <- cbind(all.svs.hmp,
                       data.frame(matrix("NN", nrow = NROW(all.svs.hmp),
                                               ncol = NCOL(hmp.after) - 11),
                                  stringsAsFactors = FALSE))
  colnames(all.svs.hmp) <- colnames(hmp.after)
  
  # filter dataframe with NN to have only SVs that were not projected for that particular cross
  all.svs.hmp.NN <- all.svs.hmp[which(!all.svs.hmp[, "pos"] %in% hmp.after[,"pos"]), ]
  all.svs.hmp.not.NN <- all.svs.hmp[which(all.svs.hmp[, "pos"] %in% hmp.after[,"pos"]), 1:2]
  
  # remove any duplicated row
  hmp.after <- hmp.after[!duplicated(hmp.after[, 1]), ]
  
  # make sure the data has the same length and same positions 
  if (all(hmp.after[, 1] == all.svs.hmp.not.NN[, 1])) {
    # merge missing data
    merged.hmp <- rbind(all.svs.hmp.NN, hmp.after)
    merged.hmp <- merged.hmp[order(merged.hmp$pos), ]
  } else {
    stop("Data have different length")
  }
  
  # write results for cross
  out.filename <- gsub("only-reseq-snps", "reseq-snps-all-crosses", filename.after.proj)
  fwrite(merged.hmp, file = out.filename, quote = FALSE, sep = "\t", na = "NA", row.names = FALSE)
  
}

cat("Done!\n\n")


### create single file in UNIX...

# cat("Creating a single file with SVs for all populations...\n")
# 
# # merge SVs only in one file
# for (i in 1:length(cross.list)) {
#   
#   cross <- cross.list[i]
#   
#   cat("  ", cross, "\n")
#   
#   filename.svs.only <- list.files(path = folder.after.proj,
#                                     pattern = "NAM_rils_projected-SVs-only",
#                                     full.names = TRUE)
#   filename.svs.only <- filename.svs.only[grep(cross, filename.svs.only)]
#   hmp.svs.only <- fread(filename.svs.only, header = TRUE, data.table = FALSE)
#   
#   if (i == 1) {
#     all.crosses.svs.only <- hmp.svs.only
#   } else {
#     all.crosses.svs.only <- cbind(all.crosses.svs.only, hmp.svs.only[12:NCOL(hmp.svs.only)])
#   }
#   
# }
# 
# # write results
# outname.final.svs <- paste0(folder.after.proj, "/NAM_rils_projected-SVs-only.all-RILs.hmp.txt")
# fwrite(all.crosses.svs.only, outname.final.svs, quote = FALSE, sep = "\t", na = NA, row.names = FALSE)
# 
# cat("Done!\n")



