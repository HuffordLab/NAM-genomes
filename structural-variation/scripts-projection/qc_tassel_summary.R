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
summary.filename <- args[1]
plot_name <- args[2]

# setwd("~/projects/sv_nams/analysis/reseq_snps_projection2")

# summary.filename = "ld/subset_high-ld-snps/tassel_summary3.txt"
# summary.filename = "ld/subset_low-ld-snps/tassel_summary3.txt"
# summary.filename = "ld/subset_random-snps/tassel_summary3.txt"
# summary.filename = "ld/tassel_summary_sv3.txt"

# plot_name <- "ld/subset_high-ld-snps/missing_snps_high.png"
# plot_name <- "ld/subset_low-ld-snps/missing_snps_low.png"
# plot_name <- "ld/subset_random-snps/missing_snps_random.png"
# plot_name <- "ld/missing_svs_subset.png"



#### libraries ----

if(!require("data.table")) install.packages("data.table")
if(!require("ggplot2")) install.packages("ggplot2")



#### plot ----

summary <- fread(summary.filename, header = TRUE, data.table = FALSE)

if (grepl("snp", summary.filename)) {
  plot_x <- "Proportion of missing per SNP"
  if (grepl("high", summary.filename)) plot_subtitle <- "Subset of SNPs in high LD to SVs"
  if (grepl("low", summary.filename)) plot_subtitle <- "Subset of SNPs in low LD to SVs"
  if (grepl("random", summary.filename)) plot_subtitle <- "Subset of random SNPs" 
} else if (grepl("sv", summary.filename)) {
  plot_x <- "Proportion of missing per SV"
  plot_subtitle <- "Subset of SVs"
}

missing_plot <- ggplot(summary, aes(x = `Proportion Missing`)) +
  geom_histogram(binwidth = 0.05) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(x = plot_x,
       y = "Count",
       fill = "Marker type",
       subtitle = plot_subtitle)

ggsave(plot = missing_plot, filename = plot_name, device = "png")
