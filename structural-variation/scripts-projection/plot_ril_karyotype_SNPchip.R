#### arguments for command line ----

args <- commandArgs(trailingOnly = TRUE)

# # help
# if (all(length(args) == 1 & args == "-h" | args == "--help")) {
#   cat("
# Description: plot karyotypes for RILs of a biparental cross
#
# Usage: ")
#   quit()
# }

# make sure the correct number of arguments are used
# you should provide 6 arguments
if (length(args) != 6) {
  stop("incorrect number of arguments provided.
       
       Usage: Rscript plot_ril_karyotypes.R [chromosomes_file] [centromeres_file] [cross] [output_folder]
       [data_filename]
       ")
}

chrs.file <- args[1]
centromeres.file <- args[2]
cross <- args[3]
output.folder <- args[4]
data.filename <- args[5]

if (grepl("--rils=", args[6])) {
  rils.list <- unlist(strsplit(args[6], split = "="))[2]
  rils.list <- unlist(strsplit(rils.list, split = ","))
} else {
  stop("Invalid list of rils. Make sure it's comma-separated")
}




# chrs.file <- "~/projects/sv_nams/analysis/qc/B73_RefGen_V4_chrm_info.txt"
# centromeres.file <- "~/projects/sv_nams/analysis/qc/centromeres_Schneider-2016-pnas_v4.bed"
# cross <- "B73xB97"
# output.folder <- "~/projects/sv_nams/analysis/qc/karyotypes/SNP_chip"
# data.filename <- "~/projects/sv_nams/data/NAM_map_and_genos-121025/hapmap/NAM_SNP_genos_raw_20090921.hmp.txt"
# rils.list <- c("Z001E0082", "Z001E0168", "Z001E0262")



#### libraries ----

if(!require("data.table")) install.packages("data.table")
if(!require("ggplot2")) install.packages("ggplot2")



#### load chromosome and centromere positions ----

# chromosomes
chrms <- read.delim(chrs.file, sep = "\t", header = TRUE)
chrms <- data.frame(chr = chrms$chr,
                    start_pos = 0,
                    end_pos = chrms$length)

# centromeres
centros <- read.delim(centromeres.file, sep = "\t", header = TRUE)
centros <- data.frame(chr = centros$chr,
                      start_pos = centros$start_pos,
                      end_pos = centros$end_pos)



#### plot karyotypes ----

# load hapmap with SNP calls for each RIL
hmp.file <- fread(data.filename, header = TRUE, data.table = FALSE)

# get parents names
parent1 <- "B73"
parent2 <- toupper(unlist(strsplit(cross, "B73x"))[2])

# get parents column numbers in gbs data
p1.col.gbs <- grep(parent1, colnames(hmp.file))
p2.col.gbs <- grep(parent2, colnames(hmp.file))

for (RIL in rils.list) {
  
  # get ril column number
  ril.col <- grep(RIL, colnames(hmp.file))
  
  # merge information of RIL of interest with respective marker positions
  geno.data <- cbind(hmp.file[, c(1,3,4)], hmp.file[, c(ril.col, p1.col.gbs, p2.col.gbs)])
  colnames(geno.data) <- c("marker", "chr", "pos", "geno", "parent1", "parent2")
  
  # transform chr column to numeric to be compatible with variables "chrms" and "centros"
  geno.data$chr <- as.numeric(geno.data$chr)
  # select missing data and non-missing data
  geno.data.not.missing <- subset(geno.data, geno != "NN")
  # get proportion of non-missing data to add in the plot
  prop.not.missing <- NROW(geno.data.not.missing) / NROW(geno.data)
  
  # check from which parent an allele came
  geno.data.not.missing.p1 <- geno.data.not.missing[which(geno.data.not.missing[, "geno"] == geno.data.not.missing[, "parent1"] &
                                                            geno.data.not.missing[, "geno"] != geno.data.not.missing[, "parent2"]), ]# &
  #geno.data.not.missing[, "parent2"] != "NN"), ]
  geno.data.not.missing.p2 <- geno.data.not.missing[which(geno.data.not.missing[, "geno"] != geno.data.not.missing[, "parent1"] &
                                                            geno.data.not.missing[, "geno"] == geno.data.not.missing[, "parent2"]), ]# &
  #geno.data.not.missing[, "parent1"] != "NN"), ]
  geno.data.not.missing.rest <- geno.data.not.missing[which(!geno.data.not.missing[, "marker"] %in% geno.data.not.missing.p1[, "marker"] &
                                                              !geno.data.not.missing[, "marker"] %in% geno.data.not.missing.p2[, "marker"]), ]
  
  geno.data.not.missing.het <- data.frame(matrix(nrow = 0, ncol = NCOL(geno.data.not.missing.rest)))
  colnames(geno.data.not.missing.het) <- colnames(geno.data.not.missing.rest)
  for (snp in 1:NROW(geno.data.not.missing.rest)) {
    # check if ril snp is a het or if it has a different allele from parents
    p1.alleles <- unlist(strsplit(geno.data.not.missing.rest[snp, "parent1"], split = ""))
    p2.alleles <- unlist(strsplit(geno.data.not.missing.rest[snp, "parent2"], split = ""))
    ril.alleles <- unlist(strsplit(geno.data.not.missing.rest[snp, "geno"], split = ""))
    if (unique(p1.alleles) %in% ril.alleles & unique(p2.alleles) %in% ril.alleles & ril.alleles[1] != ril.alleles[2]) {
      geno.data.not.missing.het <- rbind(geno.data.not.missing.het, geno.data.not.missing.rest[snp, ])
    }
  }
  
  # NROW(geno.data.not.missing.p1)+NROW(geno.data.not.missing.p2)+NROW(geno.data.not.missing.UK)+NROW(geno.data.not.missing.mono)+NROW(geno.data.not.missing.disagree)
  # NROW(geno.data.not.missing)
  
  # plot karyotype
  karyo.plot <- ggplot() +
    geom_segment(data = chrms,
                 aes(x = 0, xend = 0, y = start_pos, yend = end_pos),
                 lineend = "round", color = "Gainsboro", size = 5) +
    geom_segment(data = centros,
                 aes(x = 0, xend = 0, y = start_pos, yend = end_pos),
                 lineend = "round", color = "DimGray", size = 5) +
    geom_segment(data = geno.data.not.missing.p1,
                 aes(x = 0, xend = 0, y = pos, yend = pos + 1e6),  # increased size to be able to see the marker
                 lineend = "butt", color = "firebrick", size = 5, alpha = 0.3) +
    geom_segment(data = geno.data.not.missing.p2,
                 aes(x = 0, xend = 0, y = pos, yend = pos + 1e6),  # increased size to be able to see the marker
                 lineend = "butt", color = "#386cb0", size = 5, alpha = 0.3) +
    geom_segment(data = geno.data.not.missing.het,
                 aes(x = 0, xend = 0, y = pos, yend = pos + 1e6),  # increased size to be able to see the marker
                 lineend = "butt", color = "#fec44f", size = 5, alpha = 0.8) +
    scale_y_reverse(breaks = seq(0, 3.5e8, 0.50e8), labels = c(1, seq(50, 350, 50))) +
    scale_x_discrete(position = "top") +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.caption = element_text(size = rel(1.1), color = "DimGray"),
          axis.text = element_text(size=rel(2)),
          axis.title = element_text(size=rel(2)),
          strip.text.x = element_text(size=rel(2))) +
    facet_grid(~chr, switch = "y") +
    labs(caption = paste0(parent1, "x", parent2, " - ", gsub("RIL_", "RIL ", RIL), "\n\n",
                          "Not missing: ", round(prop.not.missing, digits = 2), "\n"),
         x = "Chromosomes", y = "Genomic positions (Mb)")
  
  # create directory for output if it doesn't exist
  if (!dir.exists(output.folder)) {
    dir.create(output.folder, recursive = TRUE)
  }
  # save plots
  karyo.name <- paste0(output.folder, "/", cross, "_", RIL, ".pdf")
  ggsave(filename = karyo.name, plot = karyo.plot, device = "pdf", width = 4.5, height = 7, units = "in")
  
}

cat("Done!\n\n")




#### credits ----

# code for karyotype plot was adapted from Carles Hernandez-Ferrer's blog:
# https://carleshf.github.io/chromosome_karyotype_plot_with_ggplot2/
