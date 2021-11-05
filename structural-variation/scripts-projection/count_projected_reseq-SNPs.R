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
       
       Usage: Rscript count_projected_SVs.R [folder_with_sv_calls] [folder_with_projected_files]
       ")
}

# assign arguments to variables
sv.folder <- args[1]
folder.after.proj <- args[2]


# sv.folder <- "~/projects/sv_nams/data/tmp"
# folder.after.proj <- "~/projects/sv_nams/analysis/reseq_snps_projection2"



#### libraries ----

if(!require("data.table")) install.packages("data.table")
if(!require("ggplot2")) install.packages("ggplot2")
if(!require("doParallel")) install.packages("doParalell")

if (detectCores() > 10) {
  num.cores <- 10
} else {
  num.cores <- detectCores()
}



#### summarize SVs of a cross ----

# create an empty df to store summary of SVs for that population
summary.projection <- data.frame(matrix(nrow = 0, ncol = 4))
colnames(summary.projection) <- c("family", "total_SVs", "avg_percent_projected", "proj_accuracy")

# create empty vector to store the percentage of SVs projected for all individuals
proj.SVs.all.rils <- c()
proj.poly.SVs.all.rils <- c()

# get list with all NAM families
cross.list <- system("ls -d ~/projects/sv_nams/data/GBS-output/tmp/B73* | xargs -n 1 basename",
                     intern = TRUE)

# get list with all SV filenames
sv.files <- list.files(sv.folder, pattern = "NAM_parents_SNPs-reseq_and_SVs-SNPs", full.names = TRUE, recursive = TRUE)
sv.files <- sv.files[grep("poly.sorted.hmp.txt", sv.files)]


# get list with all SV filenames
snp.names <- list.files(sv.folder, pattern = "all_SNP_names", full.names = TRUE)

# create empty dataframes to store percenage of projected RILs for each SNP
snp.names <- fread(snp.names, header = TRUE, data.table = FALSE)
snp.names <- snp.names[, 1]

proj.RILs.all.svs <- data.frame(matrix(NA, nrow = length(snp.names), ncol = length(cross.list) + 1), stringsAsFactors = FALSE)
colnames(proj.RILs.all.svs) <- c("sv", cross.list)
proj.RILs.all.svs$sv <- snp.names


for (cross in cross.list) {
  
  cat(cross, "\n")
  # get parents names
  parent1 <- "B73"
  parent2 <- toupper(unlist(strsplit(cross, "B73x"))[2])
  
  # get sv filename for that cross
  sv.filename <- sv.files[grep(cross, sv.files)]
  # open file with SV positions
  sv.hmp.cross <- fread(sv.filename, header = TRUE, data.table = FALSE)
  # change parent columns to upper case
  colnames(sv.hmp.cross)[12:NCOL(sv.hmp.cross)] <- toupper(colnames(sv.hmp.cross)[12:NCOL(sv.hmp.cross)])
  # rename B73 parent
  colnames(sv.hmp.cross)[grep("B73", colnames(sv.hmp.cross))] <- "B73"
  # summarize svs for that cross
  sv.hmp.cross <- sv.hmp.cross[, c("rs#", "chrom", "pos", parent1, parent2)]
  # keep only chr 1 to 10
  sv.hmp.cross <- sv.hmp.cross[which(sv.hmp.cross[, "chrom"] %in% 1:10), ]
  number.total <- NROW(sv.hmp.cross)
  
  
  
  #### summarize SV projection ----
  
  # open files after projection
  filename.after.proj <- list.files(path = folder.after.proj, pattern = paste0(cross, ".poly.projected.hmp.txt"),
                                    full.names = TRUE)
  hmp.after <- fread(filename.after.proj, header = TRUE, data.table = FALSE)
  
  
  # filter files to have only SVs
  hmp.after <- hmp.after[grep("^S[0-9]+_", hmp.after[, 1], perl = TRUE), ]
  
  # remove duplicates
  hmp.after <- hmp.after[!duplicated(hmp.after[, 1]), ]
  # make sure the same SNPs are in both files (and exclude scaffolds)
  hmp.after <- hmp.after[which(hmp.after[, 1] %in% sv.hmp.cross[, 1]), ]
  
  # count how many SVs were projected per RIL
  SVs.projected <- list()
  for (RIL in colnames(hmp.after)[12:NCOL(hmp.after)]) {
    SVs.projected[[RIL]] <- sum(hmp.after[, RIL] != "NN")
  }
  
  # count how many RILs were projected per SV
  RILs.projected <- apply(hmp.after[, 12:NCOL(hmp.after)], MARGIN = 1, function(sv) {
    return(sum(sv != "NN"))
  })
  names(RILs.projected) <- hmp.after[, 1]
  
  # create folder to store plots
  plots.folder <- paste0(folder.after.proj, "/plots_snps")
  if (!dir.exists(plots.folder)) {
    dir.create(plots.folder, recursive = TRUE)
  }
  
  # plot distribution
  proj.distribution <- as.numeric(unlist(SVs.projected))
  proj.distribution.plot <- ggplot(data.frame(dist = proj.distribution), aes(x = dist)) +
    geom_histogram() +
    labs(x = "SNPs projected",
         y = "Number of RILs")
  # add how many rils!
  ggsave(filename = paste0(plots.folder, "/", cross, "_SNPs-projection_distribution.pdf"),
         plot = proj.distribution.plot, device = "pdf")
  
  # plot distribution rils
  proj.distribution.rils <- as.numeric(RILs.projected)
  proj.distribution.rils.plot <- ggplot(data.frame(dist = proj.distribution.rils), aes(x = dist)) +
    geom_histogram() +
    labs(x = "RILs projected",
         y = "Number of SVs")
  # add how many rils!
  ggsave(filename = paste0(plots.folder, "/", cross, "_SNPs-projection_distribution_RILs.pdf"),
         plot = proj.distribution.rils.plot, device = "pdf")
  
  
  # add percent SVs projected of all individuals to vector
  proj.SVs.all.rils <- append(proj.SVs.all.rils, proj.distribution/number.total)
  
  # add percent RILs projected of all individuals to table
  if (all(proj.RILs.all.svs[which(proj.RILs.all.svs[, 1] %in% names(RILs.projected)), 1] == names(RILs.projected))) {
    
    proj.RILs.all.svs[which(proj.RILs.all.svs[, 1] %in% names(RILs.projected)), cross] <- RILs.projected / length(colnames(hmp.after)[12:NCOL(hmp.after)])
    
  } else {
    
    # keep only snps in dataset
    snps.to.keep <- names(RILs.projected)[which(names(RILs.projected) %in% proj.RILs.all.svs[, 1])]
    RILs.projected <- RILs.projected[snps.to.keep]
    sv.hmp.cross <- sv.hmp.cross[which(sv.hmp.cross[, 1] %in% snps.to.keep), ]
    hmp.after <- hmp.after[which(hmp.after[, 1] %in% snps.to.keep), ]
    
    if (all(proj.RILs.all.svs[which(proj.RILs.all.svs[, 1] %in% names(RILs.projected)), 1] == names(RILs.projected))) {
      
      # try filtering again
      proj.RILs.all.svs[which(proj.RILs.all.svs[, 1] %in% names(RILs.projected)), cross] <- RILs.projected / length(colnames(hmp.after)[12:NCOL(hmp.after)])
    
      } else {
      cat(cross, ": names don't match\n")
    }
    
  }
  
  # get accuracy
  proj.accuracy.vector <- c()
  for (chr in 1:10) {
    filename.accuracy <- list.files(path = folder.after.proj,
                                    pattern = paste0(cross, ".poly.chr-", chr, ".projected.hmp.Accuracy.txt"),
                                    full.names = TRUE)
    file.accuracy <- scan(file = filename.accuracy, what = "character", sep = "\n")
    proj.accuracy <- unlist(strsplit(file.accuracy[2], split = "\t"))
    proj.accuracy <- round(as.numeric(proj.accuracy[length(proj.accuracy)]), digits = 2)
    proj.accuracy.vector <- append(proj.accuracy.vector, proj.accuracy)
  }
  proj.accuracy.mean <- mean(proj.accuracy.vector)
  
  # get summary
  average.SVs <- round(mean(proj.distribution), digits = 1)
  percent.SVs <- round(average.SVs / number.total, digits = 2)
  
  cat(cross, ": ", average.SVs, " (", percent.SVs * 100, "%) of reseq SNPs projected on average\n", sep = "")
  
  summary.projection <- rbind(summary.projection,
                              list(family = cross,
                                   total_SVs = number.total,
                                   avg_percent_projected = percent.SVs,
                                   proj_accuracy = proj.accuracy.mean),
                              stringsAsFactors = FALSE)
  
}

# write table with projected RILs per SV
outfile.RILs.per.SV <- paste0(folder.after.proj, "/summary_projected_RILs_per_reseq-snp.txt")
fwrite(proj.RILs.all.svs, outfile.RILs.per.SV, quote = FALSE, sep = "\t", na = NA, row.names = FALSE)

# write final table
outfile <- paste0(folder.after.proj, "/summary_projection_reseq-snp.txt")
fwrite(summary.projection, outfile, quote = FALSE, sep = "\t", na = NA, row.names = FALSE)

cat("Average of ", round(mean(summary.projection$avg_percent_projected), digits = 2) * 100,
    "% reseq SNPs projected across all families (", round(mean(summary.projection$proj_accuracy) * 100, digits = 2), "% accuracy)\n\n", sep = "")

# plot final distribution of percent projected for all crosses
proj.summary.plot <- ggplot(summary.projection, aes(x = family, y = avg_percent_projected)) +
  geom_col() + 
  labs(x = "Cross",
       y = "Projected reseq SNPs (%)") +
  scale_y_continuous(labels = function(x) x*100, limits = c(0,1))+
  scale_x_discrete(labels = function(cross) gsub("B73x", "B73\nx\n", cross))

ggsave(filename = paste0(folder.after.proj, "/summary_SNPs-percent-projected_all-families.pdf"),
       plot = proj.summary.plot, device = "pdf", width = 15, height = 6)


# plot final distribution of projection accuracy for all crosses
proj.accuracy.plot <- ggplot(summary.projection, aes(x = family, y = proj_accuracy)) +
  geom_col() + 
  labs(x = "Cross",
       y = "Projection accuracy (%)") +
  scale_y_continuous(labels = function(x) x*100, limits = c(0,1))+
  scale_x_discrete(labels = function(cross) gsub("B73x", "B73\nx\n", cross))

ggsave(filename = paste0(folder.after.proj, "/summary_SNPs-projection-accuracy_all-families.pdf"),
       plot = proj.accuracy.plot, device = "pdf", width = 15, height = 6)


# plot distribution -- all rils from all families
proj.distribution.plot <- ggplot(data.frame(dist = proj.SVs.all.rils), aes(x = dist)) +
  geom_histogram() +
  labs(x = "Percent reseq SNPs projected (all RILs)",
       y = "Number of RILs") + 
  xlim(-0.1, 1.1)
# add how many rils!
ggsave(filename = paste0(folder.after.proj, "/distribution_reseq-snps-proj_all-rils.pdf"),
       plot = proj.distribution.plot, device = "pdf")


# distribution projected ril per sv
proj.RILs.all.svs.means <- rowMeans(proj.RILs.all.svs[, 2:26], na.rm = TRUE)
proj.distribution.rils.plot <- ggplot(data.frame(dist = proj.RILs.all.svs.means), aes(x = dist)) +
  geom_histogram() +
  labs(x = "Percent RILs projected",
       y = "Number of SNPs") + 
  xlim(-0.1, 1.1)
# add how many rils!
ggsave(filename = paste0(folder.after.proj, "/distribution_reseq-snps-proj_all-reseq-snps.pdf"),
       plot = proj.distribution.rils.plot, device = "pdf")

