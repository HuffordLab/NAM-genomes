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
if (length(args) != 2) {
  stop("incorrect number of arguments provided.

Usage: Rscript summary_raw_gbs.R [folder_with_raw_gbs_data] [folder_to_save_plots]
       ")
}

# assign arguments to variables
folder.raw.gbs <- args[1]
folder.for.plots <- args[2]



# folder.raw.gbs <- "~/projects/sv_nams/data/GBS-output/tmp"
# folder.for.plots <- "~/projects/sv_nams/analysis/qc/raw_gbs"



#### libraries ----

if(!require("data.table")) install.packages("data.table")
if(!require("ggplot2")) install.packages("ggplot2")




#### summarize raw gbs data ----

# create directory for plots if it doesn't exist
if (!dir.exists(folder.for.plots)) {
  dir.create(folder.for.plots, recursive = TRUE)
}

# create an empty df to store summary of SVs for that population
summary.raw.gbs <- data.frame(matrix(nrow = 0, ncol = 5))
colnames(summary.raw.gbs) <- c("family", "total_RILs", "total_SNPs", "percent_missing", "percent_polymorphic_reseq")

# create empty vector to store the average missing data per ril (all families)
percent.missing.all.rils <- c()
# create empty list to store average missing data per snp (all families)
percent.missing.all.snps <- list()


# get list with all NAM families
cross.list <- system("ls -d ~/projects/sv_nams/data/GBS-output/tmp/B73* | xargs -n 1 basename",
                     intern = TRUE)


for (cross in cross.list) {

cat("Analizing cross", cross, "...\n", sep = "")

# open files after projection
  filename <- list.files(path = paste0(folder.raw.gbs, "/", cross),
                         pattern = paste0("NAM_rils_SNPs.", cross, ".not-in-SVs.not-imputed.hmp.txt"),
                         full.names = TRUE)
  reseq.filename <- list.files(path = paste0(folder.raw.gbs, "/", cross),
                               pattern = paste0("NAM_gbs-parents_SNPs.", cross, ".not-in-SVs.reseq-overlay.hmp.txt"),
                               full.names = TRUE)
  raw.gbs <- fread(filename, header = TRUE, data.table = FALSE)
  reseq.parents <- fread(reseq.filename, header = TRUE, data.table = FALSE)
  colnames(reseq.parents)[12:NCOL(reseq.parents)]  <- toupper(colnames(reseq.parents)[12:NCOL(reseq.parents)])

  # get number of RILs and SNPs per family
  ril.cols <- grepl("Z[0-9][0-9][0-9]E", colnames(raw.gbs), perl = TRUE)
  total.RILs <- sum(ril.cols)
  total.SNPs <- NROW(raw.gbs)

  # count missing data
  ril.data <- raw.gbs[, ril.cols]
  # per ril
  missing.per.ril <- apply(X = ril.data, MARGIN = 2, FUN = function(ril) return(mean(ril == "NN")))
  percent.missing.all.rils <- append(percent.missing.all.rils, missing.per.ril)
  percent.missing.ril <- round(mean(missing.per.ril), digits = 2)
  # per snp
  missing.per.snp <- apply(X = ril.data, MARGIN = 1, FUN = function(snp) return(mean(snp == "NN")))
  percent.missing.all.snps[[cross]] <- missing.per.snp

  # plot missing data
  missing.per.ril.plot <- ggplot(data.frame(mean = missing.per.ril), aes(x = mean)) +
    geom_histogram(bins = 25) +
    labs(x = "Missing data per RIL",
         y = "Count") +
    xlim(-0.1, 1.1)
  ggsave(filename = paste0(folder.for.plots, "/", cross, "_missing-per-ril.pdf"),
         plot = missing.per.ril.plot, device = "pdf")

  missing.per.snp.plot <- ggplot(data.frame(mean = missing.per.snp), aes(x = mean)) +
    geom_histogram(bins = 25) +
    labs(x = "Missing data per SNP",
         y = "Count") +
    xlim(-0.1, 1.1)
  ggsave(filename = paste0(folder.for.plots, "/", cross, "_missing-per-snp.pdf"),
         plot = missing.per.ril.plot, device = "pdf")

  # get parents names
  parent1 <- "B73"
  parent2 <- toupper(unlist(strsplit(cross, "B73x"))[2])

  # get parents column numbers in resequencing data
  p1.col.reseq <- grep(paste0(parent1, "_"), colnames(reseq.parents))
  p2.col.reseq <- grep(paste0(parent2, "_"), colnames(reseq.parents))

  if (cross == "B73xTzi8") {
    p2.col.reseq <- grep(toupper(parent2), colnames(reseq.parents))
    # transform NA into NN
    reseq.parents[which(is.na(reseq.parents[, p2.col.reseq])), p2.col.reseq] <- "NN"
  }

  # count polymorphic markers from resequencing
  marker.type <- apply(X = reseq.parents[, c(p1.col.reseq, p2.col.reseq)],
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

  percent.poly.reseq  <- round(mean(marker.type == "poly"), digits = 2)

  # get final summary
  summary.raw.gbs <- rbind(summary.raw.gbs,
                           list(family = cross,
                                total_RILs = total.RILs,
                                total_SNPs = total.SNPs,
                                percent_missing = percent.missing.ril,
                                percent_polymorphic_reseq = percent.poly.reseq),
                              stringsAsFactors = FALSE)
  cat("Done!\n")

}

cat("Average missing data across families:", round(mean(summary.raw.gbs$percent_missing), digits = 2), "\n")

# write final table
outfile <- paste0(folder.for.plots, "/summary_raw-gbs.txt")
fwrite(summary.raw.gbs, outfile, quote = FALSE, sep = "\t", na = NA, row.names = FALSE)


# plot final distribution of missing data all crosses
proj.summary.plot <- ggplot(summary.raw.gbs, aes(x = family, y = percent_missing)) +
  geom_col() +
  labs(x = "Cross",
       y = "Average missing data (%)") +
  scale_y_continuous(labels = function(x) x*100, limits = c(0,1))+
  scale_x_discrete(labels = function(cross) gsub("B73x", "B73\nx\n", cross))

ggsave(filename = paste0(folder.for.plots, "/distribution_missing_data.pdf"),
       plot = proj.summary.plot, device = "pdf", width = 15, height = 6)


# plot distribution -- missing data from rils of all families
percent.missing.all.rils.plot <- ggplot(data.frame(dist = percent.missing.all.rils), aes(x = dist)) +
  geom_histogram() +
  labs(x = "Percent missing data per RIL (all families)",
       y = "Number of RILs") +
  xlim(-0.1, 1.1)
# add how many rils!
ggsave(filename = paste0(folder.for.plots, "/distribution_missing-per-ril_all-families.pdf"),
       plot = percent.missing.all.rils.plot, device = "pdf")


# plot distribution -- missing snps from all families
percent.missing.all.snps.df <- do.call(cbind, percent.missing.all.snps)
percent.missing.all.snps.df <- apply(percent.missing.all.snps.df, MARGIN = 1, mean)
percent.missing.all.snps.plot <- ggplot(data.frame(dist = percent.missing.all.snps.df), aes(x = dist)) +
  geom_histogram() +
  labs(x = "Percent missing data per SNP (average of all RILs)",
       y = "Number of RILs") +
  xlim(-0.1, 1.1)
# add how many rils!
ggsave(filename = paste0(folder.for.plots, "/distribution_missing-per-snp_all-families.pdf"),
       plot = percent.missing.all.snps.plot, device = "pdf")
