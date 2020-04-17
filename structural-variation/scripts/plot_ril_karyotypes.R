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
# you should provide 9 arguments
if (length(args) != 9) {
  stop("incorrect number of arguments provided.

Usage: Rscript plot_ril_karyotypes.R [chromosomes_file] [centromeres_file] [cross] [output_folder]
                                     [data_filename] [reseq_parents] [--rils] [--parents_in_data]
                                     [--overlay_reseq]
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

if (args[8] == "--parents_in_data=TRUE") {
  parents.in.data <- TRUE
} else if (args[8] == "--parents_in_data=FALSE") {
  parents.in.data <- FALSE
} else {
  stop("Invalid argument")
}

if (args[9] == "--overlay_reseq=TRUE") {
  overlay.reseq <- TRUE
} else if (args[9] == "--overlay_reseq=FALSE") {
  overlay.reseq <- FALSE
} else {
  stop("Invalid argument")
}



# # BEFORE FSFHAP - RAW GBS
# chrs.file <- "~/projects/sv_nams/analysis/qc/B73_RefGen_V4_chrm_info.txt"
# centromeres.file <- "~/projects/sv_nams/analysis/qc/centromeres_Schneider-2016-pnas_v4.bed"
# cross <- "B73xHp301"
# output.folder <- "~/projects/sv_nams/analysis/qc/karyotypes/before_imputation_raw-gbs"
# data.filename <- "B73xHp301/NAM_rils_SNPs.B73xHp301.not-in-SVs.not-imputed.hmp.txt"
# reseq.parents.filename <- "B73xHp301/NAM_gbs-parents_SNPs.B73xHp301.not-in-SVs.reseq-overlay.hmp.txt"
# random.rils <- TRUE
# parents.in.data <- TRUE
# overlay.reseq <- TRUE
# 
# # BEFORE FSFHAP - CONSENSUS
# chrs.file <- "~/projects/sv_nams/analysis/qc/B73_RefGen_V4_chrm_info.txt"
# centromeres.file <- "~/projects/sv_nams/analysis/qc/centromeres_Schneider-2016-pnas_v4.bed"
# cross <- "B73xHp301"
# output.folder <- "~/projects/sv_nams/analysis/qc/karyotypes/before_imputation_best-markers"
# data.filename <- "B73xHp301/NAM_rils_SNPs.B73xHp301.not-in-SVs.not-imputed.best-markers.hmp.txt"
# reseq.parents.filename <- "B73xHp301/NAM_gbs-parents_SNPs.B73xHp301.not-in-SVs.reseq-overlay.hmp.txt"
# #random.rils <- TRUE
# random.rils <- FALSE
# rils.list <- c("Z010E0082")
# parents.in.data <- TRUE
# overlay.reseq <- FALSE
# 
# # AFTER FSFHAP - RAW
# chrs.file <- "~/projects/sv_nams/analysis/qc/B73_RefGen_V4_chrm_info.txt"
# centromeres.file <- "~/projects/sv_nams/analysis/qc/centromeres_Schneider-2016-pnas_v4.bed"
# cross <- "B73xB97"
# output.folder <- "~/projects/sv_nams/analysis/qc/karyotypes/after_imputation_raw-gbs"
# data.filename <- "B73xB97/NAM_rils_SNPs.B73xB97.not-in-SVs.FSFHap-imputed.correct-marker-names.hmp.txt"
# reseq.parents.filename <- "B73xB97/NAM_gbs-parents_SNPs.B73xB97.not-in-SVs.reseq-overlay.hmp.txt"
# random.rils <- FALSE
# rils.list <- c("Z001E0082", "Z001E0168", "Z001E0262")
# parents.in.data <- FALSE
# overlay.reseq <- FALSE
# 
# # AFTER FSFHAP - CONSENSUS
# chrs.file <- "~/projects/sv_nams/analysis/qc/B73_RefGen_V4_chrm_info.txt"
# centromeres.file <- "~/projects/sv_nams/analysis/qc/centromeres_Schneider-2016-pnas_v4.bed"
# cross <- "B73xTzi8"
# output.folder <- "~/projects/sv_nams/analysis/qc/karyotypes/after_imputation_best-markers"
# data.filename <- "B73xTzi8/NAM_rils_SNPs.B73xTzi8.not-in-SVs.FSFHap-imputed.best-markers.hmp.txt"
# reseq.parents.filename <- "B73xTzi8/NAM_rils_SNPs.B73xTzi8.not-in-SVs.not-imputed.best-markers.hmp.txt"
# random.rils <- FALSE
# rils.list <- c("Z026E0045")
# parents.in.data <- FALSE
# overlay.reseq <- FALSE




#### libraries ----

if(!require("data.table")) install.packages("data.table")
if(!require("ggplot2")) install.packages("ggplot2")

# display warnings as they occur (if they occur)
options(warn = 1)



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

cat("\nLoading cross ", cross, "...\n", sep = "")

# load hapmap for all RILs in the RIL population
geno.data.infile <- fread(data.filename, header = TRUE, data.table = FALSE)
reseq.parents.infile <- fread(reseq.parents.filename, header = TRUE, data.table = FALSE)
# transform to all column names to upper case to avoid mismatch
colnames(geno.data.infile)[12:NCOL(geno.data.infile)]  <- toupper(colnames(geno.data.infile)[12:NCOL(geno.data.infile)])
colnames(reseq.parents.infile)[12:NCOL(reseq.parents.infile)]  <- toupper(colnames(reseq.parents.infile)[12:NCOL(reseq.parents.infile)])

cat("Done!\n")

# get parents names
parent1 <- "B73"
parent2 <- toupper(unlist(strsplit(cross, "B73x"))[2])


# if the parents are in the gbs data...
if (parents.in.data == TRUE) {
  # get parents column numbers
  p1.col.gbs <- grep(paste0(parent1, "_"), colnames(geno.data.infile))
  p2.col.gbs <- grep(paste0(parent2, "_"), colnames(geno.data.infile))
  
  # cross B73xTzi8 doesn't have gbs data on parent Tzi8. Thus, I'll have to use the resequencing data
  if (cross == "B73xTzi8") {
    p2.col.reseq <- grep(toupper(parent2), colnames(reseq.parents.infile))
    # transform NA into NN
    reseq.parents.infile[which(is.na(reseq.parents.infile[, p2.col.reseq])), p2.col.reseq] <- "NN"
    # add data to gbs file
    # simply add column if data have the same length and all positions match
    if (NROW(geno.data.infile) == NROW(reseq.parents.infile)) {
      if (all(reseq.parents.infile$pos == geno.data.infile$pos)) {
        geno.data.infile <- cbind(geno.data.infile, TZI8_PI506246 = reseq.parents.infile[, p2.col.reseq])
      } else {
        stop("Resequencing and gbs positions don't match")
      }
    } else {
      # if lengths are diferent, add only snps present in each chromosome
      geno.data.infile  <- cbind(geno.data.infile, TZI8_PI506246 = NA)
      for (chr in unique(geno.data.infile[, "chrom"])) {
        
        # subset data
        geno.data.chr <- subset(geno.data.infile, chrom == chr)
        reseq.data.chr <- subset(reseq.parents.infile, chrom == chr)
        # select positions
        SNP.positions <- geno.data.chr[, "pos"]
        # filter reseq data
        reseq.data.chr <- subset(reseq.data.chr, pos %in% SNP.positions)
        parent2.data <- reseq.data.chr[, p2.col.reseq]
        # add data to parent 2 column
        geno.data.infile[which(geno.data.infile[which(geno.data.infile[, "chrom"] == chr), "pos"] %in% SNP.positions), "TZI8_PI506246"] <- parent2.data
        
      }
    }
    
    geno.data.infile[, "TZI8_PI506246"] <- as.character(geno.data.infile[, "TZI8_PI506246"])
    
    p2.col.gbs <- grep(paste0(parent2, "_"), colnames(geno.data.infile))
    
    overlay.reseq <- FALSE
  }
  
  # make sure positions of the resequenced parents match the gbs data
  # (skip this if data comes from best markers, because parental data was already overlayed)
  if (overlay.reseq == TRUE) {
    if (all(reseq.parents.infile$pos == geno.data.infile$pos)) {
      
      cat("Getting parental genotypes from resequencing data...\n")
      # get parents column numbers in resequencing data
      p1.col.reseq <- grep(paste0(parent1, "_"), colnames(reseq.parents.infile))
      p2.col.reseq <- grep(paste0(parent2, "_"), colnames(reseq.parents.infile))
      
      # change parental data from gbs to reseq
      geno.data.infile[, p1.col.gbs] <- reseq.parents.infile[, p1.col.reseq]
      geno.data.infile[, p2.col.gbs] <- reseq.parents.infile[, p2.col.reseq]
      cat("Done!\n")
      
    } else {
      stop("Resequencing and gbs positions don't match")
    }
  }
}

# if parents are not in the gbs data
if (parents.in.data == FALSE) {
  # make sure positions of the resequenced parents match the gbs data
  if (all(reseq.parents.infile$pos == geno.data.infile$pos)) {
    
    cat("Adding parental data to file...\n")
    # get parents column numbers in resequencing data
    p1.col.reseq <- grep(paste0(parent1, "_"), colnames(reseq.parents.infile))
    p2.col.reseq <- grep(paste0(parent2, "_"), colnames(reseq.parents.infile))
    
    # change parental data from gbs to reseq
    geno.data.infile <- cbind(geno.data.infile,
                              parent1 = reseq.parents.infile[, p1.col.reseq],
                              parent2 = reseq.parents.infile[, p2.col.reseq])
    
    # get parents column numbers in gbs
    idx.p1.colname <- which(colnames(geno.data.infile) == "parent1")
    idx.p2.colname <- which(colnames(geno.data.infile) == "parent2")
    
    # correct column name of parents
    colnames(geno.data.infile)[idx.p1.colname] <- parent1
    colnames(geno.data.infile)[idx.p2.colname] <- parent2
    cat("Done!\n")
    
  } else {
    stop("Resequencing and gbs positions don't match")
  }
}


# randomly select 3 RILs per population to plot karyotype
if (random.rils == TRUE) {
  set.seed(184)
  selected.RILs <- sample(colnames(geno.data.infile[12:NCOL(geno.data.infile)]), size = 3, replace = FALSE)
  if (cross == "B73xTzi8") {
    selected.RILs <- sample(colnames(geno.data.infile[12:(NCOL(geno.data.infile)-1)]), size = 3, replace = FALSE)
  }
} else {
  selected.RILs <- rils.list
}

# plot karyotype for each selected RIL
for (RIL in selected.RILs) {
  
  cat("Plotting", RIL, "\n")
  # get ril column number
  ril.col <- grep(RIL, colnames(geno.data.infile))
  
  # merge information of RIL of interest with respective marker positions
  if (parents.in.data == TRUE) {
    geno.data <- cbind(geno.data.infile[, c(1,3,4)], geno.data.infile[, c(ril.col, p1.col.gbs, p2.col.gbs)])
  } else {
    geno.data <- cbind(geno.data.infile[, c(1,3,4)], geno.data.infile[, c(ril.col, idx.p1.colname, idx.p2.colname)])
  }
  colnames(geno.data) <- c("marker", "chr", "pos", "geno", "parent1", "parent2")
  
  # use only chromosomes from 1 to 10
  geno.data <- subset(geno.data, chr %in% 1:10)
  # transform chr column to numeric to be compatible with variables "chrms" and "centros"
  geno.data$chr <- as.numeric(geno.data$chr)
  # select missing data and non-missing data
  geno.data.not.missing <- subset(geno.data, geno != "NN")
  # get proportion of non-missing data to add in the plot
  prop.not.missing <- NROW(geno.data.not.missing) / NROW(geno.data)
  
  # check from which parent an allele came
  geno.data.not.missing.p1 <- geno.data.not.missing[which(geno.data.not.missing[, "geno"] == geno.data.not.missing[, "parent1"] &
                                                            geno.data.not.missing[, "geno"] != geno.data.not.missing[, "parent2"]), ]
  geno.data.not.missing.p2 <- geno.data.not.missing[which(geno.data.not.missing[, "geno"] != geno.data.not.missing[, "parent1"] &
                                                            geno.data.not.missing[, "geno"] == geno.data.not.missing[, "parent2"]), ]
  geno.data.not.missing.rest <- geno.data.not.missing[which(!geno.data.not.missing[, "marker"] %in% geno.data.not.missing.p1[, "marker"] &
                                                              !geno.data.not.missing[, "marker"] %in% geno.data.not.missing.p2[, "marker"]), ]
  
  # find hets
  geno.data.not.missing.het <- data.frame(matrix(nrow = 0, ncol = NCOL(geno.data.not.missing.rest)))
  colnames(geno.data.not.missing.het) <- colnames(geno.data.not.missing.rest)
  for (snp in 1:NROW(geno.data.not.missing.rest)) {
    # check if ril snp is a het or if it has a different allele from parents
    p1.alleles <- unlist(strsplit(as.character(geno.data.not.missing.rest[snp, "parent1"]), split = ""))
    p2.alleles <- unlist(strsplit(as.character(geno.data.not.missing.rest[snp, "parent2"]), split = ""))
    ril.alleles <- unlist(strsplit(as.character(geno.data.not.missing.rest[snp, "geno"]), split = ""))
    if (unique(p1.alleles) %in% ril.alleles & unique(p2.alleles) %in% ril.alleles & ril.alleles[1] != ril.alleles[2]) {
      geno.data.not.missing.het <- rbind(geno.data.not.missing.het, geno.data.not.missing.rest[snp, ])
    }
  }
  
  # plot karyotype
  karyo.plot <- ggplot() +
    geom_segment(data = chrms,
                 aes(x = 0, xend = 0, y = start_pos, yend = end_pos),
                 lineend = "round", color = "Gainsboro", size = 5) +
    geom_segment(data = centros,
                 aes(x = 0, xend = 0, y = start_pos, yend = end_pos),
                 lineend = "round", color = "DimGray", size = 5) +
    geom_segment(data = geno.data.not.missing.p1,
                 aes(x = 0, xend = 0, y = pos, yend = pos + 1e5),  # increased size to be able to see the marker
                 lineend = "butt", color = "firebrick", size = 5, alpha = 0.3) +
    geom_segment(data = geno.data.not.missing.p2,
                 aes(x = 0, xend = 0, y = pos, yend = pos + 1e5),  # increased size to be able to see the marker
                 lineend = "butt", color = "#386cb0", size = 5, alpha = 0.3) +
    geom_segment(data = geno.data.not.missing.het,
                 aes(x = 0, xend = 0, y = pos, yend = pos + 1e5),  # increased size to be able to see the marker
                 lineend = "butt", color = "#fec44f", size = 5, alpha = 0.6) +
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




#### credits ----

# code for karyotype plot was adapted from Carles Hernandez-Ferrer's blog:
# https://carleshf.github.io/chromosome_karyotype_plot_with_ggplot2/
# assign arguments to variables
