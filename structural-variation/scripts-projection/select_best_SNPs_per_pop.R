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
# you should provide 8 arguments
if (length(args) != 8) {
  stop("incorrect number of arguments provided.

Usage: Rscript select_best_SNPs_per_pop.R [cross] [data_file] [reseq_file] [folder_for_plots]
                                          [--max_missing] [--window_size] [--window_step]
                                          [--min_snps_per_window]
       ")
}

cross <- args[1]
data.filename <- args[2]
reseq.parents.filename <- args[3]
folder.for.plots <- args[4]

if (grepl("--max_missing=", args[5])) {
  threshold.missing <- as.numeric(unlist(strsplit(args[5], split = "="))[2])
} else {
  stop("Invalid argument")
}

if (grepl("--window_size=", args[6])) {
  window.size <- as.numeric(unlist(strsplit(args[6], split = "="))[2])
} else {
  stop("Invalid argument")
}

if (grepl("--window_step=", args[7])) {
  window.step <- as.numeric(unlist(strsplit(args[7], split = "="))[2])
} else {
  stop("Invalid argument")
}

if (grepl("--min_snps_per_window=", args[8])) {
  min.snps.per.window <- as.numeric(unlist(strsplit(args[8], split = "="))[2])
} else {
  stop("Invalid argument")
}

# cross <- "B73xTzi8"
# data.filename <- "B73xTzi8/NAM_rils_SNPs.B73xTzi8.not-in-SVs.not-imputed.hmp.txt"
# reseq.parents.filename <- "B73xTzi8/NAM_gbs-parents_SNPs.B73xTzi8.not-in-SVs.reseq-overlay.hmp.txt"
# folder.for.plots <- "~/projects/sv_nams/analysis/qc/filter_best_SNPs"
# threshold.missing <- 0.3
# window.size <- 15
# window.step <- 1
# min.snps.per.window <- 5




#### libraries ----

if(!require("data.table")) install.packages("data.table")
if(!require("ggplot2")) install.packages("ggplot2")
if(!require("foreach")) install.packages("foreach")
if(!require("doParallel")) install.packages("doParallel")



#### prepare datasets ----

cat("\nLoading cross ", cross, "...\n", sep = "")

# load hapmap for all RILs in the RIL population
geno.data.infile <- fread(data.filename, header = TRUE, data.table = FALSE)
reseq.parents.infile <- fread(reseq.parents.filename, header = TRUE, data.table = FALSE)

cat("Done!\n")

# transform to all column names to upper case to avoid mismatch
colnames(geno.data.infile)[12:NCOL(geno.data.infile)]  <- toupper(colnames(geno.data.infile)[12:NCOL(geno.data.infile)])
colnames(reseq.parents.infile)[12:NCOL(reseq.parents.infile)]  <- toupper(colnames(reseq.parents.infile)[12:NCOL(reseq.parents.infile)])

# get parents names
parent1 <- "B73"
parent2 <- toupper(unlist(strsplit(cross, "B73x"))[2])


# get parents column numbers in resequencing data
p1.col.reseq <- grep(paste0(parent1, "_"), colnames(reseq.parents.infile))
p2.col.reseq <- grep(paste0(parent2, "_"), colnames(reseq.parents.infile))

# cross B73xTzi8 doesn't have gbs data on parent Tzi8. Thus, I'll have to use the resequencing data
if (cross == "B73xTzi8") {
  p2.col.reseq <- grep(toupper(parent2), colnames(reseq.parents.infile))
  # transform NA into NN
  reseq.parents.infile[which(is.na(reseq.parents.infile[, p2.col.reseq])), p2.col.reseq] <- "NN"
  # add data to gbs file
  geno.data.infile <- cbind(geno.data.infile, Tzi8_PI506246 = "NN")
  colnames(geno.data.infile)[12:NCOL(geno.data.infile)]  <- toupper(colnames(geno.data.infile)[12:NCOL(geno.data.infile)])
}

# get parents column numbers in gbs data
p1.col.gbs <- grep(paste0(parent1, "_"), colnames(geno.data.infile))
p2.col.gbs <- grep(paste0(parent2, "_"), colnames(geno.data.infile))


# make sure that gbs data has the snps as reseq
# (it's an ugly filtering, but it's fast and it works)

# first create a filter with gbs positions in reseq data
gbs.filter <- c()
for (chr in 1:10) {
  chr.filter <- geno.data.infile[which(geno.data.infile[, "chrom"] == chr), "pos"] %in%
    reseq.parents.infile[which(reseq.parents.infile[, "chrom"] == chr), "pos"]
  
  gbs.filter <- append(gbs.filter, chr.filter)
}

# then, create a small hapmap file with gbs scaffolds present in reseq
# (here I'm filtering the data directly because order of scaffolds are different and some are only present in one dataset)
geno.data.scafs <- data.frame(matrix(nrow = 0, ncol = NCOL(geno.data.infile)), stringsAsFactors = FALSE)
colnames(geno.data.scafs) <- colnames(geno.data.infile)
for (scaf in unique(geno.data.infile[grep("SCAF", geno.data.infile[, "chrom"]), "chrom"])) {
  if (scaf %in% unique(reseq.parents.infile[grep("SCAF", reseq.parents.infile[, "chrom"]), "chrom"])) {
    
    gbs.scaf <- subset(geno.data.infile, chrom == scaf)
    reseq.scaf <- subset(reseq.parents.infile, chrom == scaf)
    
    gbs.scaf <- gbs.scaf[which(gbs.scaf[, "pos"] %in% reseq.scaf[, "pos"]), ]
    
    if (NROW(gbs.scaf) > 0) {
      # make sure to add only if scafold actually has shared positions between files
      geno.data.scafs <- rbind(geno.data.scafs, gbs.scaf)
    }
    
  }
}

# keep only shared positions
geno.data.infile <- geno.data.infile[gbs.filter, ]
# exclude scaffold information because above filter didn't do anything to scaffolds
geno.data.infile <- geno.data.infile[1:sum(gbs.filter), ]
# add correct scaffold info
geno.data.infile <- rbind(geno.data.infile, geno.data.scafs)

# just check if all positions match between gbs and reseq
all(geno.data.infile[, "pos"] == reseq.parents.infile[, "pos"])

# change parental data from gbs to reseq
geno.data.infile[, p1.col.gbs] <- reseq.parents.infile[, p1.col.reseq]
geno.data.infile[, p2.col.gbs] <- reseq.parents.infile[, p2.col.reseq]

# remember that this resequencing overlay file already excludes
# (i.e.change to NN) sites that disagree between gbs and reseq


cat("Checking missing data for each SNP...\n")

# for each SNP, how many RILs have it
RILs.with.SNP <- apply(X = geno.data.infile[, 15:NCOL(geno.data.infile)],
                       MARGIN = 1, FUN = function(snp) return(mean(snp != "NN")))

# create directory for plots if it doesn't exist
if (!dir.exists(folder.for.plots)) {
  dir.create(folder.for.plots, recursive = TRUE)
}

# plot distribution
RILs.with.SNP.df <- data.frame(mean = RILs.with.SNP)
RILs.with.SNP.plot <- ggplot(RILs.with.SNP.df, aes(x = mean)) +
  geom_histogram(bins = 25) +
  labs(x = "Percent of RILs with particular SNP present",
       y = "Count")
ggsave(filename = paste0(folder.for.plots, "/", cross, "_SNPs_not-missing.pdf"),
       plot = RILs.with.SNP.plot, device = "pdf")

cat("Done!\n")




#### filter by polymorphic (based on resequencing parents), and then by top not missing ----

cat("Keeping only polymorphic markers between parents...\n")

marker.type <- apply(X = geno.data.infile[, c(p1.col.gbs, p2.col.gbs)],
                     MARGIN = 1, FUN = function(snp) {

                       # get unique genotypes between parents
                       genotypes <- unique(snp)
                       genotypes <- genotypes[genotypes != "NN"]

                       if (length(genotypes) == 0) {

                         # if there is no genotype, snp is missing
                         return("missing")

                       } else if (length(genotypes) == 1) {

                         # if there is one genotype, it's monomorphic
                         # but distinguish if SNP is het
                         alleles <- unlist(strsplit(genotypes, split = ""))
                         if (alleles[1] == alleles[2]) {
                           return("mono")
                         } else {
                           return("het")
                         }

                       } else {

                         # if there are two genotypes, it's polymorphic
                         # but distiguish if one of the genotypes is het
                         p1.alleles <- unlist(strsplit(genotypes[1], split = ""))
                         p2.alleles <- unlist(strsplit(genotypes[2], split = ""))
                         if (p1.alleles[1] == p1.alleles[2] & p2.alleles[1] == p2.alleles[2]) {
                           return("poly")
                         } else {
                           return("het")
                         }

                       }
                     })

# keep only homozygous polymorphic markers between parents
geno.data.infile.poly <- geno.data.infile[which(marker.type == "poly"), ]

RILs.with.SNP.poly <- RILs.with.SNP[which(marker.type == "poly")]

cat("Done!\n")

# plot distribution
RILs.with.SNP.poly.df <- data.frame(mean = RILs.with.SNP.poly)
RILs.with.SNP.poly.plot <- ggplot(RILs.with.SNP.poly.df, aes(x = mean)) +
  geom_histogram(bins = 25) +
  labs(x = "Percent of RILs with particular SNP present",
       y = "Count")
ggsave(filename = paste0(folder.for.plots, "/", cross, "_SNPs_not-missing_polymorphic.pdf"),
       plot = RILs.with.SNP.poly.plot, device = "pdf")


cat("Keeping only SNPs that are present in more than ", threshold.missing*100, "% of RILs...\n", sep ="")

geno.data.infile.poly.top <- geno.data.infile.poly[which(RILs.with.SNP.poly > threshold.missing), ]

RILs.with.SNP.poly.top <- RILs.with.SNP.poly[which(RILs.with.SNP.poly > threshold.missing)]

cat ("Done!\n")

# plot distribution
RILs.with.SNP.poly.top.df <- data.frame(mean = RILs.with.SNP.poly.top)
RILs.with.SNP.poly.top.plot <- ggplot(RILs.with.SNP.poly.top.df, aes(x = mean)) +
  geom_histogram(bins = 25) +
  xlim(0, 1) +
  labs(x = "Percent of RILs with particular SNP present",
     y = "Count")
ggsave(filename = paste0(folder.for.plots, "/", cross, "_SNPs_not-missing_polymorphic_above-threshold.pdf"),
       plot = RILs.with.SNP.poly.top.plot, device = "pdf")


# # check if filtering was correctly done
# apply(X = geno.data.infile.poly.top[1:5, 15:NCOL(geno.data.infile.poly.top)],
#       MARGIN = 1, FUN = function(snp) return(mean(snp != "NN")))





#### sliding window approach to determine genotype call ----

cat("Getting consensus genotypes by sliding window approach...\n")

# get available cores for paralellizing
num.cores <- detectCores() - 1
# register cores for parallelizing, but limit to maximum 5 cores
if (num.cores >= 5) {
  registerDoParallel(cores = 5)
} else {
  registerDoParallel(cores = num.cores)
}


geno.data.infile.poly.top.window <- foreach(chr=1:10, .combine = rbind) %dopar% {

  # subset by chromsome
  chr.poly.window <- geno.data.infile.poly.top[which(geno.data.infile.poly.top[, "chrom"] == chr), ]
  
  # it's possible that some chromosomes don't have enough SNPs to make up a window
  # so only do sliding window if number of sites are bigger than the window size
  if (NROW(chr.poly.window) > window.size) {
    
    # for each RIL, see if call is from parent1 or parent2 (ignore missing)
    for (ril.col in 15:NCOL(chr.poly.window)) {
      
      # cat("  RIL ", colnames(chr.poly.window)[ril.col], "...", sep = "")
      
      # set up first window
      window.start <- 1
      window.stop <- window.start + (window.size - 1)
      
      # create a vector to store consenus genotype for each window
      window.consensus <- c()
      
      # use slide window approach until end of the window reaches the last SNP
      while (window.stop <= NROW(chr.poly.window)) {
        
        # get genotypes from parents and ril for that window
        window <- chr.poly.window[window.start:window.stop, c(p1.col.gbs, p2.col.gbs, ril.col)]
        
        # define from which parents the SNP in ril comes from
        window.calls <- apply(window, MARGIN = 1, FUN = function(genotypes) {
          if (genotypes[3] == "NN") {
            # if ril snp is NN
            return("missing")
          } else if (genotypes[1] == genotypes[3]) {
            # if ril snp is the same as parent1
            return("p1")
          } else if (genotypes[2] == genotypes[3]) {
            # if ril snp is the same as parent2
            return("p2")
          } else {
            # check if ril snp is a het or if it has a different allele from parents
            p1.alleles <- unlist(strsplit(genotypes[1], split = ""))
            p2.alleles <- unlist(strsplit(genotypes[2], split = ""))
            ril.alleles <- unlist(strsplit(genotypes[3], split = ""))
            if (unique(p1.alleles) %in% ril.alleles & unique(p2.alleles) %in% ril.alleles & ril.alleles[1] != ril.alleles[2]) {
              return("het")
            } else {
              return("missing")
            }
          }
        })
        
        # check if there is enough ril snps genotyped
        if (sum(window.calls != "missing") >= min.snps.per.window) {
          
          # get number of alleles for each parent (multiply by 2 because a homozygous has 2 alleles)
          n.p1.alleles <- (sum(window.calls == "p1") * 2) + (sum(window.calls == "het"))
          n.p2.alleles <- (sum(window.calls == "p2") * 2) + (sum(window.calls == "het"))
          total.alleles <- sum(window.calls != "missing") * 2
          # from those not missing, what proportion is p1, p2 and het?
          prop.p1 <- n.p1.alleles/total.alleles
          
          # define consensus based on threshold (p1: p1>0.7, p2: p1<0.3, het: 0.3<p1<0.7)
          if (prop.p1 >= 0.7) {
            # assign parent1 genotype of first snp on window to consensus
            window.consensus <- append(window.consensus, window[1, 1])
          } else if (prop.p1 <= 0.3) {
            # assign parent2 genotype of first snp on window to consensus
            window.consensus <- append(window.consensus, window[1, 2])
          } else {
            # assign het genotype to consensus
            p1.allele <- unlist(strsplit(window[1, 1], split = ""))[1]
            p2.allele <- unlist(strsplit(window[1, 2], split = ""))[1]
            window.consensus <- append(window.consensus, paste0(p1.allele, p2.allele))
          }
        } else {
          # if there is very few or none ril SNPs genotyped, consider consensus as missing data
          window.consensus <- append(window.consensus, "NN")
        }
        
        # set up the start of next window
        window.start <- window.start + window.step
        window.stop <- window.start + (window.size - 1)
      }
      
      # after that, add (window.size - 1) NNs in the consensus
      window.consensus <- append(window.consensus, rep("NN", times = (window.size - 1)))
      
      # change genotypes based on consensus
      chr.poly.window[, ril.col] <- window.consensus
      
    }
  }
  
  chr.poly.window

}
stopImplicitCluster()

cat("Done!\n")




# remove SNPs with allele frequency <0.4 or >0.6

cat("Removing SNPs with allele frequency < 0.4 or > 0.6...\n")

af.filter <- apply(geno.data.infile.poly.top.window[, 12:NCOL(geno.data.infile.poly.top.window)],
      MARGIN = 1, FUN = function(genotypes) {
        
        alleles <- paste0(genotypes, collapse = "")
        alleles <- unlist(strsplit(alleles, split = ""))
        
        allele.count <- data.frame(table(alleles))
        allele.count <- allele.count[which(allele.count[, 1] != "N"), ]
        
        allele1 <- allele.count[1, 2]
        allele2 <- allele.count[2, 2]
        
        allele.freq <- allele1 / (allele1 + allele2)
        
        return(allele.freq)
        
      })

geno.data.infile.poly.top.window.AF <- geno.data.infile.poly.top.window[which(af.filter > 0.4 & af.filter < 0.6), ]

RILs.with.SNP.poly.top.window.AF <- RILs.with.SNP.poly.top[which(af.filter > 0.4 & af.filter < 0.6)]

cat ("Done!\n")

# plot distribution
RILs.with.SNP.poly.top.window.AF.df <- data.frame(mean = RILs.with.SNP.poly.top.window.AF)
RILs.with.SNP.poly.top.window.AF.df.plot <- ggplot(RILs.with.SNP.poly.top.window.AF.df, aes(x = mean)) +
  geom_histogram(bins = 25) +
  xlim(0, 1) +
  labs(x = "Percent of RILs with particular SNP present",
       y = "Count")
ggsave(filename = paste0(folder.for.plots, "/", cross, "_SNPs_not-missing_polymorphic_above-threshold_allele-freq_filtered.pdf"),
       plot = RILs.with.SNP.poly.top.window.AF.df.plot, device = "pdf")




#### write output ----

outfile <- gsub(".hmp.txt", ".best-markers.hmp.txt", data.filename, fixed = TRUE)
fwrite(geno.data.infile.poly.top.window.AF, outfile, quote = FALSE, sep = "\t", na = NA, row.names = FALSE)
