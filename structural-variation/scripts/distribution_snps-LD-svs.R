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
if (length(args) < 2 | length(args) > 3) {
  stop("incorrect number of arguments provided.
       
       Usage: Rscript 
       ")
}

# assign arguments to variables
plink.file <- args[1]
out.dir <- args[2]
if (length(args[3]) == 0) {
  snps.to.subset <- NA
} else {
  snps.to.subset <- args[3]
}


# setwd("~/projects/sv_nams/analysis/reseq_snps_projection2")
# plink.file <- "ld/plink_results_SNPs-lowest-LD-SV_chr10.ld"
# out.dir <- "ld/subset_low-ld-snps"
# snps.to.subset <- "~/projects/sv_nams/data/subset-NAM-snps/SNPs-to-keep_subsample-low-ld_chr10.txt"

# setwd("~/projects/sv_nams/analysis/reseq_snps_projection2")
# plink.file <- "/scratch.global/della028/hirsch_lab/ld_files/NAM_rils_SNPs-reseq_and_SVs-SNPs.reseq-snps-all-crosses.chr-5.projected.duplicated-SVs-removed_ld-w-100_v2.no-tra.snp-sv.ld"
# out.dir <- "ld/subset_random-340k-snps"
# snps.to.subset <- "~/projects/sv_nams/data/subset-NAM-snps/SNPs-to-keep_subsample-340k_chr5.txt"


#### libraries ----

if(!require("data.table")) install.packages("data.table")
if(!require("ggplot2")) install.packages("ggplot2")
if(!require("doParallel")) install.packages("doParalell")


if (detectCores() > 10) {
  num.cores <- 10
} else {
  num.cores <- detectCores()
}


add.info.to.ld.df <- function(LD_results) {
  
  # add column with distance between sv and snp
  LD_results$dist_to_sv <- LD_results[, 5] - LD_results[, 2]
  
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



#### load data ----

cat("Loading data\n")

# open table with LD among markers
LD_results <- fread(plink.file, header = TRUE, data.table = FALSE)



#### analyze LD between SNPs and SVs ----

if (!dir.exists(out.dir)) dir.create(out.dir)

# get chr being analyzed
chr.plink <- LD_results[1, 1]

# get snps to keep
snps.to.keep <- fread(snps.to.subset, header = FALSE, data.table = FALSE)
snps.to.keep <- as.character(snps.to.keep[, 1])

# make sure chromosomes are the same
chr.subset <- unlist(strsplit(snps.to.keep[1], split = "_"))[1]
chr.subset <- as.numeric(gsub("S", "", chr.subset))
if (chr.subset != chr.plink) stop("Different chromosomes being used: check your input files")

# filter LD data
LD_results <- LD_results[which(LD_results[, "SNP_A"] %in% snps.to.keep | LD_results[, "SNP_B"] %in% snps.to.keep), ]


LD_results_filtered <- mclapply(snps.to.keep, FUN = function(snp, LD_results) {
  
  # subset LD results to have only the SV being parsed
  snps_notLD_with_sv <- LD_results[which(LD_results[, "SNP_A"] == snp | LD_results[, "SNP_B"] == snp), ]
  
  # select random SV that SNP is not in LD with
  if (NROW(snps_notLD_with_sv) > 0) {
    # select random SV that SNP is not in LD with
    row.to.keep <- sample(1:NROW(snps_notLD_with_sv), size = 1, replace = FALSE)
    snps_notLD_with_sv <- snps_notLD_with_sv[row.to.keep, ]
    # return random SNP and ld info
    return(snps_notLD_with_sv)
  }  
}, LD_results, mc.cores = num.cores)

# make it a data frame
LD_results_filtered <- do.call(rbind, LD_results_filtered)
# make sure they are in the right order
LD_results_filtered <- LD_results_filtered[order(LD_results_filtered$BP_A, LD_results_filtered$BP_B), ]

# remove duplicates
LD_results_filtered <- LD_results_filtered[!duplicated(LD_results_filtered[, c("SNP_A", "SNP_B")]), ]



cat("Summarizing data\n")

# write final output
outfile.name <- rev(unlist(strsplit(plink.file, "/")))[1]
outfile.name <- gsub(".ld", ".subset.ld", outfile.name, fixed = TRUE)
outfile.name <- paste0(out.dir, "/", outfile.name)
fwrite(LD_results_filtered, file = outfile.name, quote = FALSE, na = NA, sep = "\t", row.names = FALSE)


# distribution of r2 of SNPs in LD with SVs
dist.ld <- ggplot(LD_results_filtered, aes(x = R2)) +
  geom_histogram(fill = "#900721", binwidth = 0.005) +
  labs(title = paste0("LD between SVs and SNPs (chr", chr.plink, ")"),
       x = bquote("LD"~(r^2)),
       y = "Count") +
  coord_cartesian(xlim = c(0, 1)) +
  theme(title = element_text(size = 15),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15))

ggsave(dist.ld, filename = paste0(out.dir, "/dist-LD_SNPs-SVs_chr", chr.plink, ".png"), device = "png")

# add more info if missing
if (NCOL(LD_results_filtered) == 7) {
  LD_results_filtered <- add.info.to.ld.df(LD_results_filtered)
}


# summarize results
LD_results_filtered$ld_quarter <- factor(LD_results_filtered$ld_quarter,
                                         levels = c("0_to_0.25", "0.25_to_0.5", "0.5_to_0.75", "0.75_to_1"))

summary_ld <- aggregate(formula = cbind(sv_size, dist_to_sv) ~ ld_quarter, data = LD_results_filtered,
                        FUN = function(x) c(count = NROW(x), mean = mean(x), median = median(x)))
summary_ld <- do.call(data.frame, summary_ld)
summary_ld <- summary_ld[, c(1:4,6:7)]
colnames(summary_ld) <- c("ld_quarter", "sv_count", "mean_sv_size", "median_sv_size", "mean_dist_to_sv", "median_dist_to_sv")

# add counts per sv size range and ld quarter
size_range_summary <- aggregate(LD_results_filtered$sv_size,
                                by = list(LD_results_filtered$ld_quarter, LD_results_filtered$size_range),
                                FUN = function(x) NROW(x))
size_range_summary <- reshape(size_range_summary, v.names = "x", idvar = "Group.1", timevar = "Group.2", direction = "wide")
colnames(size_range_summary) <- gsub("x.", "sv_", colnames(size_range_summary), fixed = TRUE)

# write summary table
summary_ld <- cbind(summary_ld, size_range_summary[, -1])
summary_out <- paste0(out.dir, "/summary_ld_snps-svs_chr", chr.plink, ".txt")
fwrite(summary_ld, file = summary_out, quote = FALSE, na = NA, sep = "\t", row.names = FALSE)


# write final snps selected
snps.kept <- apply(LD_results_filtered, MARGIN = 1, function(row) {
  marker1 <- row["SNP_A"]
  marker2 <- row["SNP_B"]
  if (grepl(paste0("^S", chr.plink, "_"), marker1)) {
    return(marker1)
  } else {
    return(marker2)
  }
})
snps.kept <- snps.kept[!duplicated(snps.kept)]
snps.kept <- data.frame(snps.kept, stringsAsFactors = FALSE)

outfile.snps.kept <- paste0(out.dir, "/SNPs-kept_chr", chr.plink, ".txt")
fwrite(snps.kept, outfile.snps.kept, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, na = NA)

