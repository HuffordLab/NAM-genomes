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
       
       Usage: Rscript 
       ")
}

# assign arguments to variables
plink.file <- args[1]
sv.names <- args[2]
snp.names <- args[3]
out.dir.ld <- args[4]
out.dir.names <- args[5]
chr <- args[6]

# setwd("~/projects/sv_nams/analysis/reseq_snps_projection2")
# plink.file <- "/scratch.global/della028/hirsch_lab/ld_files/NAM_rils_SNPs-reseq_and_SVs-SNPs.reseq-snps-all-crosses.chr-1.projected.duplicated-SVs-removed_ld-w-100_v2.no-tra.snp-sv.ld"
# sv.names <- "~/projects/sv_nams/data/subset-NAM-snps/SVs-to-keep.missing-filter.no-tra.txt"
# snp.names <- "~/projects/sv_nams/data/subset-NAM-snps/SNPs-to-keep.missing-filter.txt"
# out.dir.ld <- "ld"
# out.dir.names <- "~/projects/sv_nams/data/subset-NAM-snps"
# chr <- 10

#### libraries ----

if(!require("data.table")) install.packages("data.table")
if(!require("doParallel")) install.packages("doParalell")


if (detectCores() > 10) {
  num.cores <- 10
} else {
  num.cores <- detectCores()
}

#### subsample ----


add.info.to.ld.df <- function(LD_results) {
  
  # add sv size
  LD_results$sv_size <- apply(LD_results, MARGIN = 1, FUN = function(info) {
    sv.index <- grep("^S[0-9]+_", info[c("SNP_A", "SNP_B")], perl = TRUE, invert = TRUE)
    sv <- info[c("SNP_A", "SNP_B")][sv.index]
    sv <- unlist(strsplit(sv, split = ".", fixed = TRUE))
    sv.size <- as.numeric(sv[4]) - as.numeric(sv[3])
    return(sv.size)
  })
  
  # add R2 quarter
  LD_results$ld_quarter <- apply(LD_results, MARGIN = 1, FUN = function(info) {
    ld <- as.numeric(info["R2"])
    if (ld >= 0 & ld <= 0.25) return("0_to_0.25")
    if (ld > 0.25 & ld <= 0.5) return("0.25_to_0.5")
    if (ld > 0.5 & ld <= 0.75) return("0.5_to_0.75")
    if (ld > 0.75 & ld <= 1.1) return("0.75_to_1")
  })
  
  # add range of sv sizes
  LD_results$size_range <- apply(LD_results, MARGIN = 1, FUN = function(info) {
    if (all(!is.na(info))) {
      sv_size <- as.numeric(info["sv_size"])
      if (sv_size >= 0 & sv_size <= 10000) return("<10kb")
      if (sv_size > 10000 & sv_size <= 100000) return("10kb-100kb")
      if (sv_size > 100000 & sv_size <= 1000000) return("100kb-1Mb")
      if (sv_size > 1000000) return(">1Mb")
    } else {
      return("NA")
    }
  })
  
  return(LD_results)
}


SVs <- fread(sv.names, header = FALSE, data.table = FALSE)
SVs <- as.character(SVs[, 1])

SNPs <- fread(snp.names, header = FALSE, data.table = FALSE)
SNPs <- as.character(SNPs[, 1])


# load one chr at a time
plink.file.chr <- gsub("chr-[0-9]+", paste0("chr-", chr), plink.file, perl = TRUE)
# open table with LD among markers
LD_results <- fread(plink.file.chr, header = TRUE, data.table = FALSE)

# subset SVs based on chr (recall that the file shouldn't have translocations)
SVs.chr <- SVs[grep(paste0(".chr", chr, "."), SVs, fixed = TRUE)]
SNPs.chr <- SNPs[grep(paste0("^S", chr, "_"), SNPs, perl = TRUE)]

#   #### qc
#   SVs.chr.notLD <- SVs.chr[which(!SVs.chr %in% LD_results[, "SNP_A"] & !SVs.chr %in% LD_results[, "SNP_B"])]
#   length(SVs.chr.notLD)
#   sv.hmp <- fread("NAM_rils_projected-SVs-only.all-RILs.duplicated-SVs-removed.chr-10.hmp.txt", header = TRUE, data.table = FALSE)
#   sv.hmp <- sv.hmp[which(sv.hmp[, 1] %in% SVs.chr.notLD), ]
#   check.mono <- apply(sv.hmp, 1, function(x) {
#     genos <- unique(x[12:length(x)])
#     genos <- unique(genos)
#     genos <- genos[genos != "NN"]
#     if (length(genos) == 0) {
#       type <- "missing"
#     } else if (length(genos) == 1) {
#       type <- "mono"
#     } else {
#       type <- "poly"
#     }
#     return(c(x[1], type))
#   })
#   check.mono <- t(check.mono)
#   colnames(check.mono) <- c("marker", "type")
#   View(check.mono[which(check.mono[, 2] == "poly"), ])
#   as.character(check.mono[which(check.mono[, 2] != "mono"), 1])
# 
#   # keep only SV-SNP LDs
#   keep_SV_SNPs_LD_only <- which(!LD_results[, 3] %in% SVs.chr & LD_results[, 6] %in% SVs.chr |
#                                   LD_results[, 3] %in% SVs.chr & !LD_results[, 6] %in% SVs.chr)
#   LD_results <- LD_results[keep_SV_SNPs_LD_only, ]
#   
#   # remove translocations
#   translocations <- apply(LD_results, MARGIN = 1, function(row) {
#     return(grepl("tra.chr", row["SNP_A"], fixed = TRUE) | grepl("tra.chr", row["SNP_B"], fixed = TRUE))
#     })
#   LD_results <- LD_results[!translocations, ]
#
# # remove duplicates
# LD_results <- LD_results[!duplicated(LD_results[, c("SNP_A", "SNP_B")]), ]

# keep only SNPs that have low % of missing data
LD_results <- LD_results[which(LD_results[, "SNP_A"] %in% SNPs.chr | LD_results[, "SNP_B"] %in% SNPs.chr), ]

# add column with distance between sv and snp
LD_results$dist_to_sv <- LD_results[, 5] - LD_results[, 2]


cat("subsetting only SNPs with highest LD to SV\n")


# create empty dataset
LD_results_highest <- data.frame(matrix(nrow = 0, ncol = NCOL(LD_results)), stringsAsFactors = FALSE)
colnames(LD_results_highest) <- colnames(LD_results)

# create vector to make sure the same snp doesn't get picked twice
snps.already.in.ld <- c()

# get closest (highest LD) snps
for (sv in SVs.chr) {
  
  # subset LD results to have only the SV being parsed
  snps_LD_with_sv <- LD_results[which(LD_results[, "SNP_A"] == sv | LD_results[, "SNP_B"] == sv), ]
  
  if (length(snps.already.in.ld) > 0) {
    snps_LD_with_sv <- snps_LD_with_sv[which(!snps_LD_with_sv[, "SNP_A"] %in% snps.already.in.ld & !snps_LD_with_sv[, "SNP_B"] %in% snps.already.in.ld), ]
  }
  
  
  if (NROW(snps_LD_with_sv) > 0) {
    
    # select only SNP with highest LD with that SV
    snps_LD_with_sv <- snps_LD_with_sv[which(snps_LD_with_sv[, "R2"] == max(snps_LD_with_sv[, "R2"])), ]
    # if there are more than one SNP with the same R2, get the closest one to the SV
    if (NROW(snps_LD_with_sv) > 1) {
      snps_LD_with_sv <- snps_LD_with_sv[which(snps_LD_with_sv[, "dist_to_sv"] == min(snps_LD_with_sv[, "dist_to_sv"])), ]
    }
    # get closest snp
    snp.selected <- apply(snps_LD_with_sv, MARGIN = 1, function(row) {
      marker1 <- row["SNP_A"]
      marker2 <- row["SNP_B"]
      if (grepl(paste0("^S", chr, "_"), marker1)) {
        return(marker1)
      } else {
        return(marker2)
      }
    })
    snps.already.in.ld <- append(snps.already.in.ld, as.character(snp.selected))
    # add closest SNP in LD with SV into new df
    LD_results_highest <- rbind(LD_results_highest, snps_LD_with_sv)
    
  }
}


#   LD_results_highest <- mclapply(SVs.chr, FUN = function(sv, LD_results) {
#     
#     # subset LD results to have only the SV being parsed
#     snps_LD_with_sv <- LD_results[which(LD_results[, "SNP_A"] == sv | LD_results[, "SNP_B"] == sv ), ]
#     
#     if (NROW(snps_LD_with_sv) > 0) {
#       
#       # select only SNP with highest LD with that SV
#       snps_LD_with_sv <- snps_LD_with_sv[which(snps_LD_with_sv[, "R2"] == max(snps_LD_with_sv[, "R2"])), ]
#       # if there are more than one SNP with the same R2, get the closest one to the SV
#       if (NROW(snps_LD_with_sv) > 1) {
#         snps_LD_with_sv <- snps_LD_with_sv[which(snps_LD_with_sv[, "dist_to_sv"] == min(snps_LD_with_sv[, "dist_to_sv"])), ]
#       }
#       
#       # add closest SNP in LD with SV into new df
#       return(snps_LD_with_sv)
#       
#     }
#     
#   }, LD_results, mc.cores = num.cores)
#   # make it a data frame
#   LD_results_highest <- do.call(rbind, LD_results_highest)
#   

# remove duplicates
LD_results_highest <- LD_results_highest[!duplicated(LD_results_highest[, c("SNP_A", "SNP_B")]), ]

# add info
LD_results_highest <- add.info.to.ld.df(LD_results_highest)

# write filtered ld table
outfile.highest <- paste0(out.dir.ld, "/plink_results_SNPs-highest-LD-SV_chr", chr, ".missing-filter.ld")
fwrite(LD_results_highest, outfile.highest, sep = "\t", quote = FALSE, row.names = FALSE, na = NA)


# get name of SNPs in highest LD to a SV
snps.highest <- apply(LD_results_highest, MARGIN = 1, function(row) {
  marker1 <- row["SNP_A"]
  marker2 <- row["SNP_B"]
  if (grepl(paste0("^S", chr, "_"), marker1)) {
    return(marker1)
  } else {
    return(marker2)
  }
})

subsample.highest <- data.frame(snps.highest, stringsAsFactors = FALSE)

outfile.subsample.highest <- paste0(out.dir.names, "/SNPs-to-keep_subsample-high-ld_chr", chr, ".missing-filter.txt")
fwrite(subsample.highest, outfile.subsample.highest, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, na = NA)



cat("subsetting only SNPs not in LD to SV\n")


#   snps.in.ld <- apply(LD_results, MARGIN = 1, function(row) {
#     
#     marker1 <- row["SNP_A"]
#     marker2 <- row["SNP_B"]
#     r2 <- as.numeric(row["R2"])
#     
#     if (grepl(paste0("^S", chr, "_"), marker1)) {
#       snp <- marker1
#     } else {
#       snp <- marker2
#     }
#     
#     if (r2 >= 0.2) return(snp)
#     
#   })

# first get names of SNPs that are in LD with an SV (R2>0.2)
snps.in.ld <- mclapply(1:NROW(LD_results), function(row, LD_results) {
  
  marker1 <- LD_results[row, "SNP_A"]
  marker2 <- LD_results[row, "SNP_B"]
  r2 <- as.numeric(LD_results[row, "R2"])
  
  if (grepl(paste0("^S", chr, "_"), marker1)) {
    snp <- marker1
  } else {
    snp <- marker2
  }
  
  if (r2 >= 0.2) return(snp)
  
}, LD_results, mc.cores = num.cores)

snps.in.ld.vector <- do.call(c, snps.in.ld)
snps.in.ld.vector <- snps.in.ld.vector[!duplicated(snps.in.ld.vector)]

# exclude such SNPs
snps.to.exclude <- which(!LD_results[, "SNP_A"] %in% snps.in.ld.vector & !LD_results[, "SNP_B"] %in% snps.in.ld.vector)
LD_results_lowest <- LD_results[snps.to.exclude, ]

# remove duplicates
LD_results_lowest <- LD_results_lowest[!duplicated(LD_results_lowest[, c("SNP_A", "SNP_B")]), ]

# add info
LD_results_lowest <- add.info.to.ld.df(LD_results_lowest)

# write filtered ld table
outfile.lowest <- paste0(out.dir.ld, "/plink_results_SNPs-lowest-LD-SV_chr", chr, ".missing-filter.ld")
fwrite(LD_results_lowest, outfile.lowest, sep = "\t", quote = FALSE, row.names = FALSE, na = NA)

# get name of snps
snps.lowest <- apply(LD_results_lowest, MARGIN = 1, function(row) {
  marker1 <- row["SNP_A"]
  marker2 <- row["SNP_B"]
  if (grepl(paste0("^S", chr, "_"), marker1)) {
    return(marker1)
  } else {
    return(marker2)
  }
})
snps.lowest <- snps.lowest[!duplicated(snps.lowest)]

# randomly select SNPs with R2<0.2
set.seed(9198)
subsample.lowest <- sample(snps.lowest, size = length(SVs.chr), replace = FALSE)
subsample.lowest <- data.frame(subsample.lowest, stringsAsFactors = FALSE)

outfile.subsample.lowest <- paste0(out.dir.names, "/SNPs-to-keep_subsample-low-ld_chr", chr, ".missing-filter.txt")
fwrite(subsample.lowest, outfile.subsample.lowest, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, na = NA)





