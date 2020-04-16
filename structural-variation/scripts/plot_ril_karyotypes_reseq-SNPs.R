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
# you should provide 7 arguments
if (length(args) != 7) {
  stop("incorrect number of arguments provided.

Usage: Rscript plot_ril_karyotypes_SVs.R [chromosomes_file] [centromeres_file] [cross] [output_folder]
                                         [data_filename] [sv_reseq_file] [--rils] [--expected-SVs]
       ")
}

chrs.file <- args[1]
centromeres.file <- args[2]
cross <- args[3]
output.folder <- args[4]
data.filename <- args[5]
reseq.parents.filename <- args[6]

if (args[7] == "--rils=random") {
  random.rils <- TRUE
} else if (grepl("--rils=", args[7])) {
  random.rils <- FALSE
  rils.list <- unlist(strsplit(args[7], split = "="))[2]
  rils.list <- unlist(strsplit(rils.list, split = ","))
} else {
  stop("Invalid list of rils. Make sure it's comma-separated")
}


# chrs.file <- "~/projects/sv_nams/analysis/qc/B73_RefGen_V4_chrm_info.txt"
# centromeres.file <- "~/projects/sv_nams/analysis/qc/centromeres_Schneider-2016-pnas_v4.bed"
# cross <- "B73xB97"
# data.filename <- "~/projects/sv_nams/analysis/reseq_snps_projection2/NAM_rils_SNPs-reseq_and_SVs-SNPs.B73xB97.poly.projected.hmp.txt"
# output.folder <- "~/projects/sv_nams/analysis/qc/karyotypes/reseq_snps_projection2"
# # random.rils  <- TRUE
# random.rils <- FALSE
# rils.list <- c("Z001E0082", "Z001E0168", "Z001E0262")
# reseq.parents.filename <- "~/projects/sv_nams/data/tmp/B73xB97/NAM_parents-reseq_SNPs.B73xB97.poly.not-in-SVs.hmp.txt"



#### libraries ----

if(!require("data.table")) install.packages("data.table")
if(!require("ggplot2")) install.packages("ggplot2")



#### functions ----

RemoveSVsOutofBounds <- function(hmp, chr_info) {
  
  filtered.hmp <- data.frame(matrix(nrow = 0, ncol = NCOL(hmp)), stringsAsFactors = FALSE)
  colnames(filtered.hmp) <- colnames(hmp)
  
  for (chr in 1:10) {
    end.pos <- chr_info[which(chr_info[, "chr"] == chr), "end_pos"]
    hmp.chr <- subset(hmp, chrom == chr)
    hmp.chr <- hmp.chr[which(hmp.chr[, "pos"] <= end.pos), ]
    filtered.hmp <- rbind(filtered.hmp, hmp.chr)
  }
  
  return(filtered.hmp)
}



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

# get names of RIL populations
cat("Plotting ", cross, "...\n", sep = "")

# load hapmap for all RILs in the RIL population
geno.data.infile <- fread(data.filename, header = TRUE, data.table = FALSE)
reseq.parents.infile <- fread(reseq.parents.filename, header = TRUE, data.table = FALSE)

# remove SVs that are outside boundaries of chromosomes (e.g. big ins/tra at the end of chr)
geno.data.infile <- RemoveSVsOutofBounds(geno.data.infile, chrms)
reseq.parents.infile <- RemoveSVsOutofBounds(reseq.parents.infile, chrms)

# transform to all column names to upper case to avoid mismatch
colnames(geno.data.infile)[12:NCOL(geno.data.infile)]  <- toupper(colnames(geno.data.infile)[12:NCOL(geno.data.infile)])
colnames(reseq.parents.infile)[12:NCOL(reseq.parents.infile)]  <- toupper(colnames(reseq.parents.infile)[12:NCOL(reseq.parents.infile)])

# get parents names
parent1 <- "B73"
parent2 <- toupper(unlist(strsplit(cross, "B73x"))[2])

# get parents column numbers in resequencing data
p1.col.reseq <- grep(paste0(parent1), colnames(reseq.parents.infile), fixed = TRUE)
p2.col.reseq <- grep(paste0(parent2), colnames(reseq.parents.infile), fixed = TRUE)


if (random.rils == TRUE) {
  # randomly select 3 RILs per population to plot karyotype
  set.seed(184)
  selected.RILs <- sample(colnames(geno.data.infile[12:NCOL(geno.data.infile)]), size = 3, replace = FALSE)
} else {
  selected.RILs <- rils.list
}

for (RIL in selected.RILs[1:2]) {
  
  # get ril column number
  ril.col <- grep(RIL, colnames(geno.data.infile))
  
  cat("  RIL", RIL, "\n")
  
  # merge information of RIL of interest with respective marker positions
  geno.data <- cbind(geno.data.infile[, c(1,3,4)], geno.data.infile[, ril.col])
  colnames(geno.data) <- c("marker", "chr", "pos", "geno")
  # use only chromosomes from 1 to 10
  geno.data <- subset(geno.data, chr %in% 1:10)
  # transform chr column to numeric to be compatible with variables "chrms" and "centros"
  geno.data$chr <- as.numeric(geno.data$chr)
  # select only reseq snps
  geno.data <- geno.data[grep("^S[0-9]+_", geno.data$marker, perl = TRUE), ]
  # remove duplicates
  geno.data <- geno.data[!duplicated(geno.data$marker), ]
  
  # add parental info
  reseq.parents <- subset(reseq.parents.infile, chrom %in% 1:10)
  # make sure they have same markers
  geno.data <- geno.data[which(geno.data[, 1] %in% reseq.parents[, 1]), ]
  geno.data <- cbind(geno.data, reseq.parents[, c(p1.col.reseq, p2.col.reseq)])
  colnames(geno.data) <- c("marker", "chr", "pos", "geno", "parent1", "parent2")
  
  # select only projected SVs
  geno.data.proj <- subset(geno.data, geno != "NN")
  geno.data.not.proj <- subset(geno.data, geno == "NN")
  
  # check from which parent an allele came
  geno.data.proj.p1 <- geno.data.proj[which(geno.data.proj[, "geno"] == geno.data.proj[, "parent1"] &
                                              geno.data.proj[, "geno"] != geno.data.proj[, "parent2"] &
                                              geno.data.proj[, "parent2"] != "NN"), ]
  geno.data.proj.p2 <- geno.data.proj[which(geno.data.proj[, "geno"] != geno.data.proj[, "parent1"] &
                                              geno.data.proj[, "geno"] == geno.data.proj[, "parent2"] &
                                              geno.data.proj[, "parent1"] != "NN"), ]
  
  geno.data.proj.het <- geno.data.proj[which(geno.data.proj[, "geno"] == "AT" | geno.data.proj[, "geno"] == "TA"), ]
  
  if (NROW(geno.data.proj.p1) == 0) {
    geno.data.proj.p1  <- rbind(geno.data.proj.p1, list(marker = NA, chr = NA, pos = NA, geno = NA, parent1 = NA, parent2 = NA))
  }
  if (NROW(geno.data.proj.p2) == 0) {
    geno.data.proj.p2  <- rbind(geno.data.proj.p2, list(marker = NA, chr = NA, pos = NA, geno = NA, parent1 = NA, parent2 = NA))
  }
  if (NROW(geno.data.proj.het) == 0) {
    geno.data.proj.het  <- rbind(geno.data.proj.het, list(marker = NA, chr = NA, pos = NA, geno = NA, parent1 = NA, parent2 = NA))
  }
  
  # plot karyotype
  karyo.plot <- ggplot() +
    geom_segment(data = chrms,
                 aes(x = 0, xend = 0, y = start_pos, yend = end_pos),
                 lineend = "round", color = "Gainsboro", size = 5) +
    geom_segment(data = centros,
                 aes(x = 0, xend = 0, y = start_pos, yend = end_pos),
                 lineend = "round", color = "DimGray", size = 5) +
    geom_segment(data = geno.data.proj.p1,
                 aes(x = 0, xend = 0, y = pos, yend = pos + 1e5),  # increased size to be able to see the marker
                 lineend = "butt", color = "#008837", size = 5, alpha = 0.3) +
    geom_segment(data = geno.data.proj.p2,
                 aes(x = 0, xend = 0, y = pos, yend = pos + 1e5),  # increased size to be able to see the marker
                 lineend = "butt", color = "#7b3294", size = 5, alpha = 0.3) +
    geom_segment(data = geno.data.proj.het,
                 aes(x = 0, xend = 0, y = pos, yend = pos + 1e5),  # increased size to be able to see the marker
                 lineend = "butt", color = "#fee08b", size = 5, alpha = 0.5) +
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
    labs(caption = paste0(cross, " - ", gsub("RIL_", "RIL ", RIL), "\n\n",
                          "Projected: ", round(NROW(geno.data.proj)/NROW(geno.data) , digits = 2), "\n"),
         x = "Chromosomes", y = "Genomic positions (Mb)")
  
  # create directory for output if it doesn't exist
  if (!dir.exists(output.folder)) {
    dir.create(output.folder, recursive = TRUE)
  }
  # save plots
  karyo.name <- paste0(output.folder, "/", cross, "_", RIL,".pdf")
  ggsave(filename = karyo.name, plot = karyo.plot, device = "pdf", width = 4.5, height = 7, units = "in")
  
}

cat("Done!\n\n")




#### credits ----

# code for karyotype plot was adapted from Carles Hernandez-Ferrer's blog:
# https://carleshf.github.io/chromosome_karyotype_plot_with_ggplot2/
# assign arguments to variables