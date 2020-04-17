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
       
       Usage: Rscript select_best_SNPs_per_pop.R [cross] [data_file] [reseq_file] [--window_size]
       [--window_step] [--min_snps_per_window]
       ")
}
cross <- args[1]
data.filename <- args[2]
reseq.parents.filename <- args[3]

if (grepl("--window_size=", args[4])) {
  window.size <- as.numeric(unlist(strsplit(args[4], split = "="))[2])
} else {
  stop("Invalid argument 4")
}

if (grepl("--window_step=", args[5])) {
  window.step <- as.numeric(unlist(strsplit(args[5], split = "="))[2])
} else {
  stop("Invalid argument 5")
}

if (grepl("--min_snps_per_window=", args[6])) {
  min.snps.per.window <- as.numeric(unlist(strsplit(args[6], split = "="))[2])
} else {
  stop("Invalid argument 6")
}

# cross <- "B73xB97"
# data.filename <- "~/projects/sv_nams/analysis/reseq_snps_projection2/NAM_rils_SNPs-reseq_and_SVs-SNPs.B73xB97.poly.chr-1.projected.hmp.txt"
# reseq.parents.filename <- "~/projects/sv_nams/analysis/reseq_snps_projection2/NAM_parents_SNPs-reseq_and_SVs-SNPs.B73xB97.poly.chr-1.sorted.hmp.txt"
# window.size <- 45
# window.step <- 1
# min.snps.per.window <- 15




#### libraries ----

if(!require("data.table")) install.packages("data.table")
if(!require("foreach")) install.packages("foreach")
if(!require("doParallel")) install.packages("doParallel")



#### sliding window approach ----

cat("\nLoading cross ", cross, "...\n", sep = "")


# get available cores for paralellizing
num.cores <- detectCores()
# register cores for parallelizing, but limit to maximum 5 cores
registerDoParallel(cores = num.cores)



# load data
geno.data.infile <- fread(data.filename, header = TRUE, data.table = FALSE)
reseq.parents.infile <- fread(reseq.parents.filename, header = TRUE, data.table = FALSE)

# get parents names
parent1 <- "B73"
parent2 <- toupper(unlist(strsplit(cross, "B73x"))[2])


# get parents column numbers in resequencing data
p1.col.reseq <- grep(parent1, colnames(reseq.parents.infile))
p2.col.reseq <- grep(parent2, colnames(reseq.parents.infile))


# make sure that gbs data has the snps as reseq
if (all(geno.data.infile[,1] != reseq.parents.infile[,1])) stop("Data have different length")


geno.data.infile.window <- foreach(ril.col=15:NCOL(geno.data.infile), .combine = cbind) %dopar% {
  
  # cat("  RIL ", colnames(geno.data.infile)[ril.col], "...", sep = "")
  
  # set up first window
  window.start <- 1
  window.stop <- window.start + (window.size - 1)
  
  # create a vector to store consenus genotype for each window
  window.consensus <- c()
  
  # use slide window approach until end of the window reaches the last SNP
  while (window.stop <= NROW(geno.data.infile)) {
    
    # get genotypes from parents and ril for that window
    window <- cbind(reseq.parents.infile[window.start:window.stop, p1.col.reseq],
                    reseq.parents.infile[window.start:window.stop, p2.col.reseq],
                    geno.data.infile[window.start:window.stop, ril.col])
    
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
  
  window.consensus
  
}
stopImplicitCluster()

# correct column names
colnames(geno.data.infile.window) <- colnames(geno.data.infile)[15:NCOL(geno.data.infile)]

# create hapmap again
geno.data.outfile <- cbind(geno.data.infile[, 1:11], geno.data.infile.window, stringsAsFactors = FALSE)

# write file
outfile <- gsub(".projected.hmp.txt", ".projected.sliding-window.hmp.txt", data.filename, fixed = TRUE)
fwrite(geno.data.outfile, outfile, quote = FALSE, sep = "\t", na = "NA", row.names = FALSE)




