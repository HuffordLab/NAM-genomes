#### arguments for command line ----

args <- commandArgs(trailingOnly = TRUE)

# # help
# if (all(length(args) == 1 & args == "-h" | args == "--help")) {
#   cat("
# Description: 
# 
# Usage: ")
#   quit()
# }

# make sure the correct number of arguments are used
# you should provide 2 arguments
if (length(args) != 1) {
  stop("incorrect number of arguments provided.
  
Usage: Rscript collapse_bionano_SVs.R [hmp_file]
       ")
}

# assign arguments to variables
hmp.file <- args[1]

# hmp.file <- "data/NAM_founders_SVs.bionano.hmp.txt"


#### libraries ----

if(!require("data.table")) install.packages("data.table")



#### collapse overlapping SVs ----

# load file
hmp <- fread(hmp.file, header = TRUE, data.table = FALSE)

# collapse svs in more than one iteration, to make sure all svs are collapsed
svs.collapsed.in.iteration <- NROW(hmp)
# count iterations
i = 1

# start loop
while(svs.collapsed.in.iteration > 0) {
  
  cat("\n---- Iteration ", i, " ----\n", sep = "")
  
  # create an empty data frame to store collapsed svs
  hmp.collapsed <- data.frame(matrix(nrow = 0, ncol = NCOL(hmp)), stringsAsFactors = FALSE)
  colnames(hmp.collapsed) <- colnames(hmp)
  # and another to keep track of those not collapsed
  hmp.not.collapsed <- data.frame(matrix(nrow = 0, ncol = NCOL(hmp)), stringsAsFactors = FALSE)
  colnames(hmp.not.collapsed) <- colnames(hmp)
  
  # analyze each chromosome separately
  for (chr in unique(hmp[, "chrom"])) {
    
    cat("Collapsing SVs from chromosome ", chr, "\n", sep = "")
    
    hmp.chr <- subset(hmp, chrom == chr)
    
    # get different types of sv in that chromosome
    sv.types <- sapply(hmp.chr[, 1], function(sv) unlist(strsplit(sv, split = ".", fixed = TRUE))[1])
    sv.types <- unique(as.character(sv.types))
    # analyze each sv type separately
    for (type in sv.types) {
      
      cat("  ", type, "\n", sep = "")
      hmp.chr.type <- subset(hmp.chr, grepl(type, hmp.chr[, 1]))
      
      # extract sv information
      sv.info <- lapply(hmp.chr.type[, 1], function(sv) unlist(strsplit(sv, split = ".", fixed = TRUE)))
      sv.info <- data.frame(do.call(rbind, sv.info), stringsAsFactors = FALSE)
      colnames(sv.info) <- c("sv_type", "chr", "start", "end")
      # make sure svs are ordered by position in chromosome
      sv.info <- sv.info[order(as.numeric(sv.info[, "start"]), as.numeric(sv.info[, "end"])), ]
      
      # create list to store overlapping svs
      overlap.list <- list()
      # loop through sv info 
      row <- 1
      while(row < NROW(sv.info)) {
        
        # get range of sv
        sv.start <- as.numeric(sv.info[row, "start"])
        sv.end <- as.numeric(sv.info[row, "end"])
        # check which svs overlap
        overlap.sv <- which(as.numeric(sv.info[, "start"]) <= sv.end)
        # make sure to exclude previous svs
        overlap.sv <- overlap.sv[!overlap.sv <= row]
        
        if (length(overlap.sv) > 0) {
          # get names of overlapping svs
          overlap.sv.names <- apply(sv.info[seq(row, max(overlap.sv), by = 1), ], MARGIN = 1, function(sv) {
            return(paste0(sv, collapse = "."))
          })
          # now add to list
          overlap.list[[length(overlap.list) + 1]] <- overlap.sv.names
          # go to row corresponding to the last overlapping sv
          row <- max(overlap.sv) + 1
        } else {
          # or go to next line if no duplicates
          row <- row + 1
        }
      }
      
      # create vector with non-overlapping svs
      non.overlapping.svs <- hmp.chr.type[which(!hmp.chr.type[, 1] %in% as.character(unlist(overlap.list))), 1]
      # add non-overlaping svs to final data frame
      if (length(non.overlapping.svs) > 0) {
        hmp.not.collapsed <- rbind(hmp.not.collapsed, hmp.chr.type[which(hmp.chr.type[, 1] %in% non.overlapping.svs), ],
                                   stringsAsFactors = FALSE)
      }
      
      # collapse calls for each group of overlapping svs
      if (length(overlap.list) > 0) {
        
        for (sv.group in 1:length(overlap.list)) {
          # get all rows that have that position
          hmp.chr.type.overlap <- hmp.chr.type[which(hmp.chr.type[, 1] %in% overlap.list[[sv.group]]), ]
          # collapse genotype calls for all lines
          collapsed.SVs <- sapply(X = hmp.chr.type.overlap[, 12:NCOL(hmp.chr.type.overlap)],
                                  FUN = function(geno) {
                                    # if there's only one unique genotype called in all duplicates
                                    if (length(unique(geno)) == 1) {
                                      # return that genotype
                                      return(unique(geno))
                                      # if there's more than one genotype call
                                    } else {
                                      # given that bionano calls are only AA or TT, if there's more
                                      # than one unique genotype, it means that there is a SV there
                                      return("TT")
                                    }
                                  })
          
          # define boundaries of collapsed sv based on the smallest start and largest end positions
          collapsed.sv.pos <- sapply(overlap.list[[sv.group]], function(sv) {
            return(unlist(strsplit(sv, split = ".", fixed = TRUE))[3:4])
          })
          collapsed.sv.start <- min(as.numeric(collapsed.sv.pos[1, ]))
          collapsed.sv.end <- max(as.numeric(collapsed.sv.pos[2, ]))
          # rewrite name of sv
          collapsed.sv.name <- paste0(type, ".chr", chr, ".", collapsed.sv.start, ".", collapsed.sv.end)
          # get middle position of sv
          collapsed.middle.pos <- ceiling((collapsed.sv.start + collapsed.sv.end) / 2)
          # define hmp structure
          hmp.info.sv <- data.frame(collapsed.sv.name, "A/T", chr, collapsed.middle.pos,
                                    hmp.chr.type.overlap[1, 5:11])
          # create a single row with collapsed sv information and collapsed calls
          hmp.chr.type.overlap.collapsed <- cbind(hmp.info.sv, as.list(collapsed.SVs), stringsAsFactors = FALSE)
          colnames(hmp.chr.type.overlap.collapsed) <- colnames(hmp.chr.type.overlap)
          
          # add this row to the final data frame
          hmp.collapsed <- rbind(hmp.collapsed, hmp.chr.type.overlap.collapsed, stringsAsFactors = FALSE)
          
        }
      }
    }
  }
  
  # merge hmps of collapsed and not collapsed svs
  hmp.final <- rbind(hmp.collapsed, hmp.not.collapsed, stringsAsFactors = FALSE)
  # sort by chr and position
  hmp.final <- hmp.final[order(hmp.final[, "chrom"], hmp.final[, "pos"]), ]
  hmp.final[, 1] <- as.character(hmp.final[, 1])
  
  # check how many svs were collapsed at this iteration
  svs.collapsed.in.iteration <- NROW(hmp) - NROW(hmp.final)
  cat("Number of collapsed SVs in this iteration: ", svs.collapsed.in.iteration, "\n", sep = "")
  # keep track of iteration
  i = i + 1
  # make final hmp as the initial hmp for next iteration
  hmp <- hmp.final
  
}

# write final table
outfile <- gsub("hmp.txt", "collapsed.hmp.txt", hmp.file)
fwrite(hmp.final, file = outfile, sep = "\t", quote = FALSE, row.names = FALSE, na = NA)
