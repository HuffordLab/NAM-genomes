if(!require("data.table")) install.packages("data.table")
if(!require("foreach")) install.packages("foreach")
if(!require("doParallel")) install.packages("doParallel")


# sv.folder <- "~/projects/sv_nams/data"
folder.after.proj <- "~/projects/sv_nams/analysis/projection"


# get list with all NAM families
cross.list <- system("ls -d ~/projects/sv_nams/data/GBS-output/tmp/B73* | xargs -n 1 basename",
                     intern = TRUE)

# get a list with all sv positions per chromosome for all crosses
all.proj.sv.pos <- list()

cat("Getting list with SV positions for all populations...\n")

for (cross in cross.list) {

  # get parents names
  parent1 <- "B73"
  parent2 <- toupper(unlist(strsplit(cross, "B73x"))[2])

  # load hapmap after projection
  filename.after.proj <- list.files(path = folder.after.proj,
                                    pattern = "projected.hmp.txt",
                                    full.names = TRUE)
  filename.after.proj <- filename.after.proj[grep(cross, filename.after.proj)]
  hmp.after <- fread(filename.after.proj, header = TRUE, data.table = FALSE)

  # filter hmp by svs and select only chromosome and positions
  sv.info <- hmp.after[grep(".", hmp.after[, 1], fixed = TRUE), c(1, 3, 4)]

  for (chr in unique(sv.info[, "chrom"])) {

    sv.info.chr <- subset(sv.info, chrom == chr)

    if (!as.character(chr) %in% names(all.proj.sv.pos)) {
      all.proj.sv.pos[[as.character(chr)]] <- sv.info.chr[, "pos"]
    } else {
      all.proj.sv.pos[[as.character(chr)]] <- append(all.proj.sv.pos[[as.character(chr)]], sv.info.chr[, "pos"])
      all.proj.sv.pos[[as.character(chr)]] <- sort(unique(all.proj.sv.pos[[as.character(chr)]]))
    }

  }

}

cat("Done!\n\n")

cat("Creating file with SVs only for each cross...\n")

for (cross in cross.list) {

  cat("  ", cross, "\n")

  # load all SVs called for all RILs
  all.svs.hmp <- fread("~/projects/sv_nams/data/NAM_founders_SVs.hmp.txt", header = TRUE, data.table = FALSE)
  all.svs.hmp <- all.svs.hmp[1:11]

  # load SVs projected for a cross
  filename.after.proj <- list.files(path = folder.after.proj,
                                    pattern = "projected.hmp.txt",
                                    full.names = TRUE)
  filename.after.proj <- filename.after.proj[grep(cross, filename.after.proj)]
  hmp.after <- fread(filename.after.proj, header = TRUE, data.table = FALSE)
  hmp.after.svs <- hmp.after[grep(".", hmp.after[, 1], fixed = TRUE), ]


  # get available cores for paralellizing
  num.cores <- detectCores() - 1
  # register cores for parallelizing, but limit to maximum 10 cores
  if (num.cores >= 10) {
    registerDoParallel(cores = 10)
  } else {
    registerDoParallel(cores = num.cores)
  }


  only.proj.svs <- foreach(chr=1:10, .combine = rbind) %dopar% {

    # filter all SVs by chr and by projected SVs
    only.proj.svs.chr <- subset(all.svs.hmp, chrom == chr)
    sv.filter <- which(only.proj.svs.chr[, "pos"] %in% all.proj.sv.pos[[as.character(chr)]])
    only.proj.svs.chr <- only.proj.svs.chr[sv.filter, ]

    # filter projected SVs by chr
    hmp.after.svs.chr <- subset(hmp.after.svs, chrom == chr)

    # cbind NN matrix with column number = number of genotypes in a cross
    only.proj.svs.chr <- cbind(only.proj.svs.chr,
                               data.frame(matrix("NN", nrow = NROW(only.proj.svs.chr),
                                                 ncol = NCOL(hmp.after.svs.chr) - 11),
                                          stringsAsFactors = FALSE))
    colnames(only.proj.svs.chr) <- colnames(hmp.after.svs.chr)

    # filter dataframe with NN to have only SVs that were not projected for that particular cross
    only.proj.svs.chr.NN <- only.proj.svs.chr[which(!only.proj.svs.chr[, 4] %in% hmp.after.svs.chr[,4]), ]
    only.proj.svs.chr.not.NN <- only.proj.svs.chr[which(only.proj.svs.chr[, 4] %in% hmp.after.svs.chr[,4]), ]

    for (row in 1:NROW(only.proj.svs.chr.not.NN)) {

      proj.filter <- which(hmp.after.svs.chr[, 1] == only.proj.svs.chr.not.NN[row, 1] & hmp.after.svs.chr[, 4] == only.proj.svs.chr.not.NN[row, 4])
      if (length(proj.filter) == 1) {
        only.proj.svs.chr.not.NN[row, ] <- as.character(hmp.after.svs.chr[proj.filter, ])
      } else if (length(proj.filter) > 1)  {
        only.proj.svs.chr.not.NN[row, ] <- as.character(hmp.after.svs.chr[proj.filter[1], ])
      }

    }

    # then rbind above df to table with projected SVs
    only.proj.svs.chr <- rbind(only.proj.svs.chr.NN, only.proj.svs.chr.not.NN)
    only.proj.svs.chr[, "pos"] <- as.integer(only.proj.svs.chr[, "pos"])
    # finally sort df by position
    only.proj.svs.chr <- only.proj.svs.chr[order(only.proj.svs.chr[, "pos"]), ]

    # return results from chr
    only.proj.svs.chr

  }
  stopImplicitCluster()


  # keep only RILs
  ril.columns <- grep("Z[0-9][0-9][0-9]E", colnames(only.proj.svs), perl = TRUE)
  only.proj.svs <- cbind(only.proj.svs[, 1:11], only.proj.svs[, ril.columns])

  # write results for cross
  out.filename <- paste0(folder.after.proj, "/NAM_rils_projected-SVs-only.", cross, ".hmp.txt")
  fwrite(only.proj.svs, file = out.filename, quote = FALSE, sep = "\t", na = "NA", row.names = FALSE)

}

cat("Done!\n\n")

cat("Creating a single file with SVs for all populations...\n")

# merge SVs only in one file
for (i in 1:length(cross.list)) {

  cross <- cross.list[i]

  cat("  ", cross, "\n")

  filename.svs.only <- list.files(path = folder.after.proj,
                                    pattern = "NAM_rils_projected-SVs-only",
                                    full.names = TRUE)
  filename.svs.only <- filename.svs.only[grep(cross, filename.svs.only)]
  hmp.svs.only <- fread(filename.svs.only, header = TRUE, data.table = FALSE)

  if (i == 1) {
    all.crosses.svs.only <- hmp.svs.only
  } else {
    all.crosses.svs.only <- cbind(all.crosses.svs.only, hmp.svs.only[12:NCOL(hmp.svs.only)])
  }

}

# write results
outname.final.svs <- paste0(folder.after.proj, "/NAM_rils_projected-SVs-only.all-RILs.hmp.txt")
fwrite(all.crosses.svs.only, outname.final.svs, quote = FALSE, sep = "\t", na = NA, row.names = FALSE)

cat("Done!\n")
