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
subset.folder <- args[1]
plot_name <- args[2]

# setwd("~/projects/sv_nams/analysis/reseq_snps_projection2")

# subset.folder = "ld/subset_high-ld-snps"
# subset.folder = "ld/subset_low-ld-snps"
# subset.folder = "ld/subset_random-snps"

# plot_name <- "ld/subset_high-ld-snps/dist-LD_SNPs-SVs_high.png"
# plot_name <- "ld/subset_low-ld-snps/dist-LD_SNPs-SVs_low.png"
# plot_name <- "ld/subset_random-snps/dist-LD_SNPs-SVs_random.png"




#### libraries ----

if(!require("data.table")) install.packages("data.table")
if(!require("ggplot2")) install.packages("ggplot2")



#### plot ----

ld.files <- list.files(path = subset.folder, pattern = ".subset.ld", full.names = TRUE)

ld.df <- data.frame(stringsAsFactors = FALSE)
for (chr in 1:10) {
  ld.chr <- ld.files[grep(paste0("chr", chr, "."), ld.files, fixed = TRUE)]
  if (length(ld.chr) == 0) ld.chr <- ld.files[grep(paste0("chr-", chr, "."), ld.files, fixed = TRUE)]
  ld.chr <- fread(ld.chr, header = TRUE, data.table = FALSE)
  ld.df <- rbind(ld.df, ld.chr)
}


if (grepl("high", subset.folder)) plot_subtitle <- "Subset of SNPs in high LD to SVs"
if (grepl("low", subset.folder)) plot_subtitle <- "Subset of SNPs in low LD to SVs"
if (grepl("random", subset.folder)) plot_subtitle <- "Subset of random SNPs"


# distribution of r2 of SNPs in LD with SVs
dist.ld <- ggplot(ld.df, aes(x = R2)) +
  geom_histogram(fill = "#900721", binwidth = 0.005) +
  labs(x = bquote("LD"~(r^2)),
       y = "Count",
       subtitle = plot_subtitle) +
  coord_cartesian(xlim = c(0, 1)) +
  theme(title = element_text(size = 15),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15))

ggsave(plot = dist.ld, filename = plot_name, device = "png")
