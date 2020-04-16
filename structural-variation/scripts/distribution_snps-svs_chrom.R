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
snp_with_sv_data <- args[1]
plot_name <- args[2]

# setwd("~/projects/sv_nams/analysis/reseq_snps_projection2")

# snp_with_sv_data = "ld/subset_high-ld-snps/SNPs-kept_chr1.txt"
# snp_with_sv_data = "ld/subset_low-ld-snps/SNPs-kept_chr1.txt"
# snp_with_sv_data = "ld/subset_random-snps/SNPs-kept_chr1.txt"

# plot_name <- "ld/subset_high-ld-snps/distribution_snps_chrom_high.png"
# plot_name <- "ld/subset_low-ld-snps/distribution_snps_chrom_low.png"
# plot_name <- "ld/subset_random-snps/distribution_snps_chrom_random.png"





#### libraries ----

if(!require("data.table")) install.packages("data.table")
if(!require("ggplot2")) install.packages("ggplot2")



#### plot ----

# get the correct genotypic data
geno_data <- data.frame(stringsAsFactors = FALSE)
for (chr in 1:10) {
  geno_data_chr <- gsub("_chr[0-9]+", paste0("_chr", chr), snp_with_sv_data, perl = TRUE)
  geno_data_chr <- fread(geno_data_chr, header = FALSE, data.table = FALSE)
  geno_data <- rbind(geno_data, geno_data_chr)
}
geno_data_chr <- NULL
geno_data <- as.character(geno_data[, 1])


# create list of SNPs
list_of_SNP_IDs <- geno_data[grep("^S", geno_data, perl = TRUE)]
SNP_pos <- strsplit(list_of_SNP_IDs, split = "_", fixed = TRUE)
SNP_pos <- data.frame(do.call(rbind, SNP_pos), stringsAsFactors = FALSE)
SNP_pos[, 1] <- gsub("^S", "", SNP_pos[, 1], perl = TRUE)
SNP_pos <- data.frame(marker_type = "snp", chrom = as.numeric(SNP_pos[, 1]), pos = as.numeric(SNP_pos[, 2]), marker = "SNP",
                      stringsAsFactors = FALSE)


# create list of SVs (non-translocations)
list_of_SV_IDs <- geno_data[grep("^S|^tra", geno_data, perl = TRUE, invert = TRUE)]
if (length(list_of_SV_IDs) > 0) {
  SV_pos <- strsplit(list_of_SV_IDs, split = ".", fixed = TRUE)
  SV_pos <- data.frame(do.call(rbind, SV_pos), stringsAsFactors = FALSE)
  SV_pos[, 2] <- as.numeric(gsub("chr", "", SV_pos[, 2], perl = TRUE))
  sv_middle_point <- apply(SV_pos[, 3:4], MARGIN = 1, function(sv) {
    start <- as.numeric(sv[1])
    end <- as.numeric(sv[2])
    middle <- ceiling(mean(c(start, end)))
    return(middle)
  })
  SV_pos <- data.frame(marker_type = SV_pos[, 1], chrom = SV_pos[, 2], pos = sv_middle_point, marker = "SV",
                       stringsAsFactors = FALSE)
  SV_pos <- SV_pos[order(SV_pos$chrom, SV_pos$pos), ]
  
  dist_svs <- ggplot(SV_pos, aes(x = pos, fill = marker)) + 
    geom_histogram(alpha = 0.5, position = "identity", binwidth = 5000000) +
    facet_wrap(~chrom, scales = "free_x") +
    scale_x_continuous(labels = function(x) x/1000000) +
    labs(x = "Position (Mb)",
         y = "Count",
         fill = "Marker type")
  
  # merge datasets
  markers_pos <- rbind(SNP_pos, SV_pos, stringsAsFactors = FALSE)
  
} else {
  
  markers_pos <- SNP_pos
  
}


# markers_pos <- markers_pos[which(markers_pos[, "chrom"] == 1), ]
# ggplot(markers_pos, aes(x = pos, fill = marker_type)) +
#   geom_density(alpha = 0.5) +
#   facet_wrap(~chrom, scales = "free")


if (grepl("high", snp_with_sv_data)) plot_subtitle <- "Subset of SNPs in high LD to SVs"
if (grepl("low", snp_with_sv_data)) plot_subtitle <- "Subset of SNPs in low LD to SVs"
if (grepl("random", snp_with_sv_data)) plot_subtitle <- "Subset of random SNPs"


dist_hist <- ggplot(markers_pos, aes(x = pos, fill = marker)) + 
  geom_histogram(position = "identity", binwidth = 5000000) +
  facet_wrap(~chrom, scales = "free_x") +
  scale_x_continuous(labels = function(x) x/1000000) +
  scale_color_manual(values = "gray50", aesthetics = "fill") +
  labs(x = "Position (Mb)",
       y = "Count",
       fill = "Marker type",
       subtitle = plot_subtitle)

ggsave(plot = dist_hist, filename = plot_name, device = "png")

# dist_densi <- ggplot(markers_pos, aes(x = pos, fill = marker)) + 
#   geom_density(alpha = 0.5, position = "identity", binwidth = 5000000) +
#   facet_wrap(~chrom, scales = "free_x") +
#   scale_x_continuous(labels = function(x) x/1000000) +
#   labs(x = "Position (Mb)",
#        y = "Count",
#        fill = "Marker type")
