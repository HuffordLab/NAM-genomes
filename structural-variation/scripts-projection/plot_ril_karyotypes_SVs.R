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
                                         [data_filename] [sv_reseq_file] [--rils]
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
# data.filename <- "~/projects/sv_nams/analysis/imputation/NAM_rils_SVs-SNPs.B73xB97.best-markers.projected.hmp.txt"
# output.folder <- "~/projects/sv_nams/analysis/qc/karyotypes/SV_projection/"
# # random.rils  <- TRUE
# random.rils <- FALSE
# rils.list <- c("Z001E0082", "Z001E0168", "Z001E0262")
# reseq.parents.filename <- "~/projects/sv_nams/data/NAM_founders_SVs.hmp.txt"



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

# get names of RIL populations
cat("Plotting ", cross, "...\n", sep = "")

# load hapmap for all RILs in the RIL population
geno.data.infile <- fread(data.filename, header = TRUE, data.table = FALSE)
reseq.parents.infile <- fread(reseq.parents.filename, header = TRUE, data.table = FALSE)
# transform to all column names to upper case to avoid mismatch
colnames(geno.data.infile)[12:NCOL(geno.data.infile)]  <- toupper(colnames(geno.data.infile)[12:NCOL(geno.data.infile)])
colnames(reseq.parents.infile)[12:NCOL(reseq.parents.infile)]  <- toupper(colnames(reseq.parents.infile)[12:NCOL(reseq.parents.infile)])

# get parents names
parent1 <- "B73"
parent2 <- toupper(unlist(strsplit(cross, "B73x"))[2])

# get parents column numbers in resequencing data
p1.col.reseq <- grep(paste0(parent1, "."), colnames(reseq.parents.infile), fixed = TRUE)
p2.col.reseq <- grep(paste0(parent2, "."), colnames(reseq.parents.infile), fixed = TRUE)


if (random.rils == TRUE) {
  # randomly select 3 RILs per population to plot karyotype
  set.seed(184)
  selected.RILs <- sample(colnames(geno.data.infile[12:NCOL(geno.data.infile)]), size = 3, replace = FALSE)
} else {
  selected.RILs <- rils.list
}

for (RIL in selected.RILs) {
  
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
  
  # select all SVs
  geno.data.all.SVs <- geno.data[grep("^del|^dup|^inv|^tra|^ins", geno.data$marker, perl = TRUE), ]
  # add parental info
  reseq.parents <- subset(reseq.parents.infile, chrom %in% 1:10)
  geno.data.all.SVs <- cbind(geno.data.all.SVs, reseq.parents[, c(p1.col.reseq, p2.col.reseq)])
  colnames(geno.data.all.SVs) <- c("marker", "chr", "pos", "geno", "parent1", "parent2")

  # select SVs per type
  geno.data.del <- geno.data[grep("del.", geno.data$marker, fixed = TRUE), ]
  geno.data.dup <- geno.data[grep("dup.", geno.data$marker, fixed = TRUE), ]
  geno.data.inv <- geno.data[grep("inv.", geno.data$marker, fixed = TRUE), ]
  geno.data.tra <- geno.data[grep("tra.", geno.data$marker, fixed = TRUE), ]
  geno.data.ins <- geno.data[grep("ins.", geno.data$marker, fixed = TRUE), ]
  
  # NROW(geno.data.del)+NROW(geno.data.dup)+NROW(geno.data.inv)+NROW(geno.data.ins)+NROW(geno.data.tra)
  # NROW(geno.data.all.SVs)
  
  # select only projected SVs
  geno.data.all.SVs.proj <- subset(geno.data.all.SVs, geno != "NN")
  geno.data.all.SVs.not.proj <- subset(geno.data.all.SVs, geno == "NN")
  
  # check from which parent an allele came
  geno.data.all.SVs.proj.p1 <- geno.data.all.SVs.proj[which(geno.data.all.SVs.proj[, "geno"] == geno.data.all.SVs.proj[, "parent1"] &
                                                              geno.data.all.SVs.proj[, "geno"] != geno.data.all.SVs.proj[, "parent2"] &
                                                              geno.data.all.SVs.proj[, "parent2"] != "NN"), ]
  geno.data.all.SVs.proj.p2 <- geno.data.all.SVs.proj[which(geno.data.all.SVs.proj[, "geno"] != geno.data.all.SVs.proj[, "parent1"] &
                                                              geno.data.all.SVs.proj[, "geno"] == geno.data.all.SVs.proj[, "parent2"] &
                                                              geno.data.all.SVs.proj[, "parent1"] != "NN"), ]
  
  geno.data.all.SVs.proj.het <- geno.data.all.SVs.proj[which(geno.data.all.SVs.proj[, "geno"] == "AT" | geno.data.all.SVs.proj[, "geno"] == "TA"), ]
  
  if (NROW(geno.data.all.SVs.proj.p1) == 0) {
    geno.data.all.SVs.proj.p1  <- rbind(geno.data.all.SVs.proj.p1, list(marker = NA, chr = NA, pos = NA, geno = NA, parent1 = NA, parent2 = NA))
  }
  if (NROW(geno.data.all.SVs.proj.p2) == 0) {
    geno.data.all.SVs.proj.p2  <- rbind(geno.data.all.SVs.proj.p2, list(marker = NA, chr = NA, pos = NA, geno = NA, parent1 = NA, parent2 = NA))
  }
  if (NROW(geno.data.all.SVs.proj.het) == 0) {
    geno.data.all.SVs.proj.het  <- rbind(geno.data.all.SVs.proj.het, list(marker = NA, chr = NA, pos = NA, geno = NA, parent1 = NA, parent2 = NA))
  }
  
  # select projected SVs by type
  geno.data.del.proj <- subset(geno.data.del, geno != "NN")
  geno.data.dup.proj <- subset(geno.data.dup, geno != "NN")
  geno.data.inv.proj <- subset(geno.data.inv, geno != "NN")
  geno.data.tra.proj <- subset(geno.data.tra, geno != "NN")
  geno.data.ins.proj <- subset(geno.data.ins, geno != "NN")
  
  # plot karyotype
  karyo.plot <- ggplot() +
    geom_segment(data = chrms,
                 aes(x = 0, xend = 0, y = start_pos, yend = end_pos),
                 lineend = "round", color = "Gainsboro", size = 5) +
    geom_segment(data = centros,
                 aes(x = 0, xend = 0, y = start_pos, yend = end_pos),
                 lineend = "round", color = "DimGray", size = 5) +
    geom_segment(data = geno.data.all.SVs.proj.p1,
                 aes(x = 0, xend = 0, y = pos, yend = pos + 1e5),  # increased size to be able to see the marker
                 lineend = "butt", color = "#008837", size = 5, alpha = 0.3) +
    geom_segment(data = geno.data.all.SVs.proj.p2,
                 aes(x = 0, xend = 0, y = pos, yend = pos + 1e5),  # increased size to be able to see the marker
                 lineend = "butt", color = "#7b3294", size = 5, alpha = 0.3) +
    geom_segment(data = geno.data.all.SVs.proj.het,
                 aes(x = 0, xend = 0, y = pos, yend = pos + 1e5),  # increased size to be able to see the marker
                 lineend = "butt", color = "#fee08b", size = 5, alpha = 0.5) +
#     geom_segment(data = geno.data.del,
#                  aes(x = 0, xend = 0, y = pos, yend = pos + 1e5),  # increased size to be able to see the marker
#                  lineend = "butt", color = "#7b3294", size = 5, alpha = 0.3) +
#     geom_segment(data = geno.data.dup,
#                  aes(x = 0, xend = 0, y = pos, yend = pos + 1e5),  # increased size to be able to see the marker
#                  lineend = "butt", color = "#008837", size = 5, alpha = 0.3) +
#     geom_segment(data = geno.data.tra,
#                  aes(x = 0, xend = 0, y = pos, yend = pos + 1e5),  # increased size to be able to see the marker
#                  lineend = "butt", color = "#7b3294", size = 5, alpha = 0.3) +
#     geom_segment(data = geno.data.inv,
#                  aes(x = 0, xend = 0, y = pos, yend = pos + 1e5),  # increased size to be able to see the marker
#                  lineend = "butt", color = "#008837", size = 5, alpha = 0.3) +
#     geom_segment(data = geno.data.ins,
#                  aes(x = 0, xend = 0, y = pos, yend = pos + 1e5),  # increased size to be able to see the marker
#                  lineend = "butt", color = "#008837", size = 5, alpha = 0.3) +
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
                          "Projected: ", round(NROW(geno.data.all.SVs.proj)/NROW(geno.data.all.SVs) , digits = 2), "\n"),
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