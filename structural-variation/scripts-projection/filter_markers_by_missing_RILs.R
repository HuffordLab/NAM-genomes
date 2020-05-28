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
# you should provide 3 arguments
if (length(args) != 6) {
  stop("incorrect number of arguments provided.
  
Usage: Rscript count_projected_SVs.R [folder_with_sv_calls] [folder_with_projected_files]
       ")
}

# assign arguments to variables
RILs.per.SV.filename <- args[1]
svs.no.duplicates.filename <- args[2]
RILs.per.SNP.filename <- args[3]
snps.after.missing.filter <- args[4]
threshold <- args[5]
out.folder <- args[6]


# RILs.per.SV.filename <- "~/projects/sv_nams/analysis/projection/summary_projected_RILs_per_sv.txt"
# svs.no.duplicates.filename <- "~/projects/sv_nams/analysis/projection/SV_names_after_removing_duplicates.txt"
# RILs.per.SNP.filename <- "~/projects/sv_nams/analysis/reseq_snps_projection2/summary_projected_RILs_per_reseq-snp.txt"
# snps.after.missing.filter <- "~/projects/sv_nams/analysis/reseq_snps_projection2/ld/missing_data_filter/SNPs_low_missing-data.txt"
# threshold <- 0.8
# out.folder <- "~/projects/sv_nams/analysis/reseq_snps_projection2/ld/missing_data_filter"


#### libraries ----

if(!require("data.table")) install.packages("data.table")




#### filter SVs by missing data in RILs ---- 

# load SVs
RILs.per.SV.df <- fread(RILs.per.SV.filename, header = TRUE, data.table = FALSE)
# keep only SVs with more than 80% data (for thoses with information on founder)
RILs.per.SV.means <- data.frame(sv = RILs.per.SV.df[, 1],
                                mean = rowMeans(RILs.per.SV.df[, 2:26], na.rm = TRUE),
                                stringsAsFactors = FALSE)

svs.no.duplicates <- fread(svs.no.duplicates.filename, header = TRUE, data.table = FALSE)
svs.no.duplicates <- svs.no.duplicates[ , 1]


# keep only SVs that are not duplicates with more than 80% data
SVs.to.keep <- RILs.per.SV.means[which(RILs.per.SV.means[, 2] >= threshold), 1]
SVs.to.keep <- SVs.to.keep[SVs.to.keep %in% svs.no.duplicates]
SVs.to.keep <- SVs.to.keep[!duplicated(SVs.to.keep)]
# remove translocations
SVs.to.keep <- SVs.to.keep[grep("^tra", SVs.to.keep, perl = TRUE, invert = TRUE)]

# load SNPs
RILs.per.SNP.df <- fread(RILs.per.SNP.filename, header = TRUE, data.table = FALSE)
# keep only SNPs with more than 80% data across all families
snps.after.missing.filter <- fread(snps.after.missing.filter, header = FALSE, data.table = FALSE)
snps.after.missing.filter <- snps.after.missing.filter[, 1]
RILs.per.SNP.df <- RILs.per.SNP.df[which(RILs.per.SNP.df[, 1] %in% snps.after.missing.filter), ]
# keep only SNPs with more than 80% data (for thoses with information on founder)
RILs.per.SNP.means <- data.frame(sv = RILs.per.SNP.df[, 1],
                                 mean = rowMeans(RILs.per.SNP.df[, 2:26], na.rm = TRUE),
                                 stringsAsFactors = FALSE)

SNPs.to.keep <- RILs.per.SNP.means[which(RILs.per.SNP.means[, 2] >= threshold), 1]
SNPs.to.keep <- SNPs.to.keep[!duplicated(SNPs.to.keep)]


# write markers to keep by chr
for (chr in 1:10) {
  
  SVs.to.keep.chr <- SVs.to.keep[grep(paste0(".chr", chr, "."), SVs.to.keep, fixed = TRUE)] 
  SNPs.to.keep.chr <- SNPs.to.keep[grep(paste0("^S", chr, "_"), SNPs.to.keep, perl = TRUE)]
  
  markers.to.keep <- append(SVs.to.keep.chr, SNPs.to.keep.chr)
  markers.to.keep <- data.frame(markers.to.keep)
  
  outfile <- paste0(out.folder, "/markers_to_keep_chr-", chr, "_missing_filter.txt")
  fwrite(markers.to.keep, outfile, quote = FALSE, row.names = FALSE, col.names = FALSE)
  
}


