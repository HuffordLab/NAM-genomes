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


# sv.folder <- "~/projects/sv_nams/data"
# folder.after.proj <- "~/projects/sv_nams/analysis/projection"
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



#### functions ----

SelectOnlySVs <- function(hmp.snps.svs, sv.info) {
  
  hmp.filtered <- data.frame(matrix(nrow = 0, ncol = NCOL(hmp.snps.svs)))
  colnames(hmp.filtered) <- colnames(hmp.snps.svs)
  
  for (chr in unique(sv.info[, "chrom"])) {
    
    sv.info.chr <- subset(sv.info, chrom == chr)
    positions.with.SVs <- sv.info.chr[, "pos"]
    
    hmp.snps.svs.chr <- subset(hmp.snps.svs, chrom == chr)
    hmp.snps.svs.chr.SVs <- hmp.snps.svs.chr[which(hmp.snps.svs.chr[, "pos"] %in% positions.with.SVs), ]
    hmp.filtered <- rbind(hmp.filtered, hmp.snps.svs.chr.SVs)
    
  }
  
  return(hmp.filtered)
}


#### summarize SVs of a cross ----

# create an empty df to store summary of SVs for that population
summary.projection <- data.frame(matrix(nrow = 0, ncol = 7))
colnames(summary.projection) <- c("family", "total_SVs", "missing_SVs", "polymorphic_SVs", "avg_percent_projected", "avg_percent_projected_poly", "proj_accuracy")

# get list with all NAM families
cross.list <- system("ls -d ~/projects/sv_nams/data/GBS-output/tmp/B73* | xargs -n 1 basename",
                     intern = TRUE)

# get list with all SV filenames
sv.files <- list.files(sv.folder, pattern = "NAM_founders_SVs_", full.names = TRUE)

# create empty vector to store the percentage of SVs projected for all individuals
proj.SVs.all.rils <- c()
proj.poly.SVs.all.rils <- c()

# create empty dataframes to store percenage of projected RILs for each SV
sv.names <- fread(sv.files[1], header = TRUE, data.table = FALSE)
sv.names <- sv.names[, 1]

proj.RILs.all.svs <- data.frame(matrix(NA, nrow = length(sv.names), ncol = length(cross.list) + 1), stringsAsFactors = FALSE)
colnames(proj.RILs.all.svs) <- c("sv", cross.list)
proj.RILs.all.svs$sv <- sv.names

# proj.poly.RILs.all.svs <- data.frame(matrix(NA, nrow = length(sv.names), ncol = length(cross.list) + 1), stringsAsFactors = FALSE)
# colnames(proj.poly.RILs.all.svs) <- c("sv", cross.list)
# proj.poly.RILs.all.svs$sv <- sv.names



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
  
  # get information of donors (parents of a cross) for projections
  filename.donors <- list.files(path = sv.folder, pattern = paste0("NAM_parents_SVs-SNPs.", cross, ".sorted.hmp.txt"),
                                full.names = TRUE)
  hmp.donors <- fread(filename.donors, header = TRUE, data.table = FALSE)
  # change parent columns to upper case
  colnames(hmp.donors)[12:NCOL(hmp.donors)] <- toupper(colnames(hmp.donors)[12:NCOL(hmp.donors)])
  # keep only SVs
  hmp.donors.svs <- hmp.donors[which(hmp.donors[, 1] %in% sv.hmp.cross[, 1]), ]
  
  # mising SVS in both parents
  missing.SVs.both <- hmp.donors.svs[, parent1] == "NN" & hmp.donors.svs[, parent2] == "NN"
  # missing SVs in one of the parents
  missing.SVs <- hmp.donors.svs[, parent1] == "NN" | hmp.donors.svs[, parent2] == "NN"
  
  # hmp without missing SVs
  hmp.donors.svs.not.missing <- hmp.donors.svs[!missing.SVs, ]
  # pololymorphic SVs
  hmp.donors.svs.poly <- hmp.donors.svs.not.missing[which(hmp.donors.svs.not.missing[, parent1] != hmp.donors.svs.not.missing[, parent2]), ]
  
  number.total <- NROW(hmp.donors.svs)
  number.missing <- sum(hmp.donors.svs[, parent2] == "NN")  # missing in nonB73 parent
  number.poly <- NROW(hmp.donors.svs.poly)
  
  
  #### summarize SV projection ----

  # open files after projection
  filename.after.proj <- list.files(path = folder.after.proj, pattern = paste0(cross, ".poly.projected.hmp.txt"),
                                    full.names = TRUE)
  hmp.after <- fread(filename.after.proj, header = TRUE, data.table = FALSE)

  # filter files to have only SVs
  hmp.after.filtered <- hmp.after[which(hmp.after[, 1] %in% sv.hmp.cross[, 1]), ]
  
  # count how many SVs were projected per RIL
  SVs.projected <- list()
  for (RIL in colnames(hmp.after.filtered)[12:NCOL(hmp.after.filtered)]) {
    SVs.projected[[RIL]] <- sum(hmp.after.filtered[, RIL] != "NN")
  }
  
  # count how many RILs were projected per SV
  RILs.projected <- apply(hmp.after.filtered[, 12:NCOL(hmp.after.filtered)], MARGIN = 1, function(sv) {
    return(sum(sv != "NN"))
  })
  names(RILs.projected) <- hmp.after.filtered[, 1]
  
  # filter data frame to have only projected polymorphic SVs
  hmp.after.filtered.poly <- data.frame(matrix(nrow = 0, ncol = NCOL(hmp.after.filtered)))
  colnames(hmp.after.filtered.poly) <- colnames(hmp.after.filtered)
  for (chr in unique(hmp.donors.svs.poly[, "chrom"])) {
    # subset data
    hmp.donors.svs.poly.chr <- subset(hmp.donors.svs.poly, chrom == chr)
    hmp.after.filtered.chr <- subset(hmp.after.filtered, chrom == chr)
    # select positions
    poly.SVs.ids <- hmp.donors.svs.poly.chr[, "rs#"]
    poly.SVs.positions <- hmp.donors.svs.poly.chr[, "pos"]
    hmp.after.filtered.chr <- subset(hmp.after.filtered.chr, `rs#` %in% poly.SVs.ids & pos %in% poly.SVs.positions)
    # append to new df
    hmp.after.filtered.poly <- rbind(hmp.after.filtered.poly, hmp.after.filtered.chr)
  }
  # count how many polymorphic SVs were projected per RIL
  SVs.projected.poly <- list()
  for (RIL in colnames(hmp.after.filtered.poly)[12:NCOL(hmp.after.filtered.poly)]) {
    SVs.projected.poly[[RIL]] <- sum(hmp.after.filtered.poly[, RIL] != "NN")
  }
  # count how many RILs were projected per SV
  RILs.projected.poly <- apply(hmp.after.filtered.poly[, 12:NCOL(hmp.after.filtered.poly)], MARGIN = 1, function(sv) {
    return(sum(sv != "NN"))
  })
  names(RILs.projected.poly) <- hmp.after.filtered.poly[, 1]
  
  # create folder to store plots
  plots.folder <- paste0(folder.after.proj, "/plots_svs")
  if (!dir.exists(plots.folder)) {
    dir.create(plots.folder, recursive = TRUE)
  }
  
  # plot distribution
  proj.distribution <- as.numeric(unlist(SVs.projected))
  proj.distribution.plot <- ggplot(data.frame(dist = proj.distribution), aes(x = dist)) +
    geom_histogram() +
    labs(x = "SVs projected",
         y = "Number of RILs")
    # add how many rils!
  ggsave(filename = paste0(plots.folder, "/", cross, "_SVs-projection_distribution.pdf"),
         plot = proj.distribution.plot, device = "pdf")
  
  # plot distribution polymorphic
  proj.distribution.poly <- as.numeric(unlist(SVs.projected.poly))
  proj.distribution.poly.plot <- ggplot(data.frame(dist = proj.distribution.poly), aes(x = dist)) +
    geom_histogram() +
    labs(x = "Polymorphic SVs projected",
         y = "Number of RILs")
  # add how many rils!
  ggsave(filename = paste0(plots.folder, "/", cross, "_SVs-projection_distribution_poly.pdf"),
         plot = proj.distribution.poly.plot, device = "pdf")
  
  # plot distribution rils
  proj.distribution.rils <- as.numeric(RILs.projected)
  proj.distribution.rils.plot <- ggplot(data.frame(dist = proj.distribution.rils), aes(x = dist)) +
    geom_histogram() +
    labs(x = "RILs projected",
         y = "Number of SVs")
  # add how many rils!
  ggsave(filename = paste0(plots.folder, "/", cross, "_SVs-projection_distribution_RILs.pdf"),
         plot = proj.distribution.rils.plot, device = "pdf")
  
  # plot distribution polymorphic rils
  proj.distribution.poly.rils <- as.numeric(RILs.projected.poly)
  proj.distribution.poly.rils.plot <- ggplot(data.frame(dist = proj.distribution.poly.rils), aes(x = dist)) +
    geom_histogram() +
    labs(x = "RILs projected",
         y = "Number of polymorphic SVs")
  # add how many rils!
  ggsave(filename = paste0(plots.folder, "/", cross, "_SVs-projection_distribution_poly-RILs.pdf"),
         plot = proj.distribution.poly.rils.plot, device = "pdf")
  
  
  # add percent SVs projected of all individuals to vector
  proj.SVs.all.rils <- append(proj.SVs.all.rils, proj.distribution/NROW(hmp.after.filtered))
  proj.poly.SVs.all.rils <- append(proj.poly.SVs.all.rils, proj.distribution.poly/NROW(hmp.after.filtered.poly))
  
  # proj.RILs.all.svs <- append(proj.RILs.all.svs, proj.distribution.rils/length(colnames(hmp.after.filtered)[12:NCOL(hmp.after.filtered)]))
  # proj.poly.RILs.all.svs <- append(proj.poly.RILs.all.svs, proj.distribution.poly.rils/length(colnames(hmp.after.filtered)[12:NCOL(hmp.after.filtered)]))
  
  #   proj.RILs.all.svs[1:10000, cross] <- sapply(proj.RILs.all.svs[1:10000, "sv"], function(sv) {
  #     
  #     if (sv %in% names(RILs.projected)) {
  #       return(RILs.projected[[sv]])
  #     } else {
  #       return(NA)
  #     }
  #     
  #   })
  
  proj.RILs.all.svs.cross <- mclapply(proj.RILs.all.svs[, "sv"], function(sv) {
    
    if (sv %in% names(RILs.projected)) {
      return(RILs.projected[[sv]])
    } else {
      return(NA)
    }
    
  }, mc.cores = num.cores)
  proj.RILs.all.svs.cross <- do.call(c, proj.RILs.all.svs.cross)
  
  # modify dataframe
  proj.RILs.all.svs[, cross] <- proj.RILs.all.svs.cross / length(colnames(hmp.after.filtered)[12:NCOL(hmp.after.filtered)])
  
  
  # for (sv in 1:length(RILs.projected.poly)) {
  #  perc.ril.proj.poly <- as.numeric(RILs.projected.poly[sv]) / length(colnames(hmp.after.filtered)[12:NCOL(hmp.after.filtered)])
  #  proj.poly.RILs.all.svs[[names(RILs.projected.poly)[sv]]] <- append(proj.poly.RILs.all.svs[[names(RILs.projected.poly)[sv]]], perc.ril.proj.poly)
  # }
  
  
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
  percent.SVs <- round(average.SVs / NROW(hmp.after.filtered), digits = 2)
  average.SVs.poly <- round(mean(proj.distribution.poly), digits = 1)
  percent.SVs.poly <- round(average.SVs.poly / NROW(hmp.after.filtered.poly), digits = 2)
  
  
  cat(cross, ": ", average.SVs, " (", percent.SVs * 100, "%) of SVs projected on average\n", sep = "")
  
  summary.projection <- rbind(summary.projection,
                              list(family = cross,
                                   total_SVs = number.total,
                                   missing_SVs = number.missing,
                                   polymorphic_SVs = number.poly,
                                   avg_percent_projected = percent.SVs,
                                   avg_percent_projected_poly = percent.SVs.poly,
                                   proj_accuracy = proj.accuracy.mean),
                              stringsAsFactors = FALSE)
  
}


# write table with projected RILs per SV
outfile.RILs.per.SV <- paste0(folder.after.proj, "/summary_projected_RILs_per_sv.txt")
fwrite(proj.RILs.all.svs, outfile.RILs.per.SV, quote = FALSE, sep = "\t", na = NA, row.names = FALSE)


# write final table
outfile <- paste0(folder.after.proj, "/summary_projection_sv.txt")
fwrite(summary.projection, outfile, quote = FALSE, sep = "\t", na = NA, row.names = FALSE)

cat("Average of ", round(mean(summary.projection$avg_percent_projected), digits = 2) * 100,
    "% SVs projected across all families (", round(mean(summary.projection$proj_accuracy), digits = 2), "% accuracy)\n\n", sep = "")

# plot final distribution of percent projected for all crosses
proj.summary.plot <- ggplot(summary.projection, aes(x = family, y = avg_percent_projected)) +
  geom_col() + 
  labs(x = "Cross",
       y = "Projected SVs (%)") +
  scale_y_continuous(labels = function(x) x*100, limits = c(0,1))+
  scale_x_discrete(labels = function(cross) gsub("B73x", "B73\nx\n", cross))

ggsave(filename = paste0(folder.after.proj, "/summary_SVs-percent-projected_all-families.pdf"),
       plot = proj.summary.plot, device = "pdf", width = 15, height = 6)


# plot final distribution of percent projected for all crosses
proj.summary.plot.poly <- ggplot(summary.projection, aes(x = family, y = avg_percent_projected_poly)) +
  geom_col() + 
  labs(x = "Cross",
       y = "Projected polymorphic SVs (%)") +
  scale_y_continuous(labels = function(x) x*100, limits = c(0,1))+
  scale_x_discrete(labels = function(cross) gsub("B73x", "B73\nx\n", cross))

ggsave(filename = paste0(folder.after.proj, "/summary_SVs-percent-poly-projected_all-families.pdf"),
       plot = proj.summary.plot.poly, device = "pdf", width = 15, height = 6)


# plot final distribution of projection accuracy for all crosses
proj.accuracy.plot <- ggplot(summary.projection, aes(x = family, y = proj_accuracy)) +
  geom_col() + 
  labs(x = "Cross",
       y = "Projection accuracy (%)") +
  scale_y_continuous(labels = function(x) x*100, limits = c(0,1))+
  scale_x_discrete(labels = function(cross) gsub("B73x", "B73\nx\n", cross))

ggsave(filename = paste0(folder.after.proj, "/summary_SVs-projection-accuracy_all-families.pdf"),
       plot = proj.accuracy.plot, device = "pdf", width = 15, height = 6)


# plot distribution -- all rils from all families
proj.distribution.plot <- ggplot(data.frame(dist = proj.SVs.all.rils), aes(x = dist)) +
  geom_histogram() +
  labs(x = "Percent SVs projected (all RILs)",
       y = "Number of RILs") + 
  xlim(-0.1, 1.1)
# add how many rils!
ggsave(filename = paste0(folder.after.proj, "/distribution_SVs-proj_all-rils.pdf"),
       plot = proj.distribution.plot, device = "pdf")


# plot distribution polymorphic -- all rils from all families
proj.distribution.poly.plot <- ggplot(data.frame(dist = proj.poly.SVs.all.rils), aes(x = dist)) +
  geom_histogram() +
  labs(x = "Percent polymorphic SVs projected (All RILs)",
       y = "Number of RILs") +
  xlim(-0.1, 1.1)
# add how many rils!
ggsave(filename = paste0(folder.after.proj, "/distribution_SVs-poly-proj_all-rils.pdf"),
       plot = proj.distribution.poly.plot, device = "pdf")


# distribution projected ril per sv
proj.RILs.all.svs.means <- rowMeans(proj.RILs.all.svs[, 2:26], na.rm = TRUE)
proj.distribution.rils.plot <- ggplot(data.frame(dist = proj.RILs.all.svs.means), aes(x = dist)) +
  geom_histogram() +
  labs(x = "Percent RILs projected",
       y = "Number of SVs") + 
  xlim(-0.1, 1.1)
# add how many rils!
ggsave(filename = paste0(folder.after.proj, "/distribution_SVs-proj_all-SVs.pdf"),
       plot = proj.distribution.rils.plot, device = "pdf")



# proj.distribution.poly.rils.plot <- ggplot(data.frame(dist = proj.poly.RILs.all.svs), aes(x = dist)) +
#   geom_histogram() +
#   labs(x = "Percent RILs projected",
#        y = "Number of polymorphic SVs") + 
#   xlim(-0.1, 1.1)
# # add how many rils!
# ggsave(filename = paste0(folder.after.proj, "/distribution_poly-sv-proj_all-SVs.pdf"),
#        plot = proj.distribution.poly.rils.plot, device = "pdf")
# 


