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
if (length(args) != 2) {
  stop("incorrect number of arguments provided.
       
       Usage: Rscript remove_duplicated_SVs.R [parents_file] [rils_file]
       ")
}

# assign arguments to variables
sv.parents.file <- args[1]
sv.rils.file <- args[2]

# sv.parents.file <- "data/NAM_founders_SVs.hmp.txt"
# sv.rils.file <- "analysis/projection/NAM_rils_projected-SVs-only.all-RILs.final.hmp.txt"


#### libraries ----

if(!require("data.table")) install.packages("data.table")
if(!require("doParallel")) install.packages("doParallel")




#### functions ----

SortMixedColumn <- function(hmp) {
  
  # subset by chomosomes with numbers
  hmp.num <- hmp[which(hmp[, "chrom"] %in% 1:10), ]
  hmp.num[, "chrom"] <- as.numeric(hmp.num[, "chrom"])
  # subset by chromosomes with characters
  hmp.chr <- hmp[which(!hmp[, "chrom"] %in% 1:10), ]
  
  # order each subset separately
  hmp.num <- hmp.num[order(hmp.num[, "chrom"], hmp.num[, "pos"]), ]
  hmp.chr <- hmp.chr[order(hmp.chr[, "chrom"], hmp.chr[, "pos"]), ]
  
  # merge the two subsets
  hmp.sorted <- rbind(hmp.num, hmp.chr)
  
  return(hmp.sorted)
  
}




#### load data ----

cat("Loading parental data\n", sep = "")
sv.parents <- fread(sv.parents.file, header = TRUE, data.table = FALSE)




#### remove duplicated SVs in parental data ----

# set maximum number of cores to be used in parallel
if (detectCores() > 10) {
  num.cores <- 10
} else {
  num.cores <- detectCores()
}

# remove duplicated SVs from parental data (except translocations)
sv.parents.tra <- sv.parents[grepl("^tra", sv.parents[,1], perl = TRUE), ]
sv.parents.no.tra <- sv.parents[!grepl("^tra", sv.parents[,1], perl = TRUE), ]
# create df to store results
sv.parents.collapsed.no.tra <- data.frame(matrix(nrow = 0, ncol = NCOL(sv.parents)))
colnames(sv.parents.collapsed.no.tra) <- colnames(sv.parents)

cat("Removing any duplicated SV from parental data\n", sep = "")
# find and collapse duplicates per chromosome
for (chr in unique(sv.parents.no.tra[, "chrom"])) {
  
  cat("  chromosome ", chr, "\n", sep = "")
  
  # subset by chromosome
  sv.parents.no.tra.chr <- sv.parents.no.tra[which(sv.parents.no.tra[, "chrom"] == chr), ]
  duplicates <- sv.parents.no.tra.chr[, 1][duplicated(sv.parents.no.tra.chr[, 1])]
  # since some rows can be repeated more than 2 times, need to get unique duplicated values
  duplicates <- unique(duplicates)
  
  # get collapsed rows
  collapsed.rows <- mclapply(duplicates, FUN = function(dup.id, chr.df) {
    
    # get all rows that have that position
    chr.df.dup <- chr.df[which(chr.df[, 1] == dup.id), ]
    
    # collapse genotype calls for all lines
    collapsed.SVs <- sapply(X = chr.df.dup[, 12:NCOL(chr.df.dup)],
                             FUN = function(geno) {
                               # if there's only one unique genotype called in all duplicates
                               if (length(unique(geno)) == 1) {
                                 # return that genotype
                                 return(unique(geno))
                                 
                                 # if there's more than one genotype call
                               } else {
                                 # check how many other non-missing genotypes there are
                                 not.NN.geno <- geno[which(geno != "NN")]
                                 # if there's only one unique genotype other than NN
                                 if (length(unique(not.NN.geno)) == 1) {
                                   # return that genotype
                                   return(unique(not.NN.geno))
                                 } else {
                                   # if there's more than one, it means there's a disagreement between calls,
                                   # so return "NN"
                                   return("NN")
                                 }
                               }
                             })
    
    # save information about SNPs
    SV.info <- chr.df.dup[1, 1:11]
    
    filtered.row <- cbind(SV.info, as.list(collapsed.SVs), stringsAsFactors = FALSE)
    
    # return row of dataframe
    return(filtered.row)
    
  }, sv.parents.no.tra.chr, mc.cores = num.cores)
  
  # transform lists into dataframe
  collapsed.df <- do.call(rbind, collapsed.rows)
  # remove all duplicated rows from main df
  sv.parents.no.tra.chr.no.dups <- sv.parents.no.tra.chr[which(!sv.parents.no.tra.chr[, 1] %in% duplicates), ]
  # merge collapsed.df to main df without duplicates
  sv.parents.no.tra.chr.filtered <- rbind(sv.parents.no.tra.chr.no.dups, collapsed.df)
  # sort dataframeby chr and position
  sv.parents.no.tra.chr.filtered <- sv.parents.no.tra.chr.filtered[order(sv.parents.no.tra.chr.filtered[, "chrom"],
                                                     sv.parents.no.tra.chr.filtered[, "pos"]), ]
  
  # append to df with other chromosomes
  sv.parents.collapsed.no.tra <- rbind(sv.parents.collapsed.no.tra, sv.parents.no.tra.chr.filtered)
  
}

# add translocations back
sv.parents <- rbind(sv.parents.collapsed.no.tra, sv.parents.tra)

# since chrom column has both numbers and characters, order them separately and then merge
sv.parents <- SortMixedColumn(sv.parents)




#### find and remove overlapping SVs ----

# create a list with any value just to start the loop
# (this value will be replaced inside the loop)
svs.to.remove <- c("start")
# start counting the number of iterations
i = 1


while (length(svs.to.remove) > 0) {
  
  cat("\n--- Iteration #", i, " ---\n", sep = "")
  
  # create a df with all sv information
  sv.info.df <- apply(sv.parents[, c(1,3,4)], MARGIN = 1, function(sv) {
    
    # split information from sv name
    info <- unlist(strsplit(sv[1], split = ".", fixed = TRUE))
    type <- info[1]
    if (type == "tra") {
      # use coordinates of where the translocation starts in B73
      chr <- paste0("chr", sv[2])
      start <- as.numeric(sv[3])
      end <- as.numeric(sv[3])
    } else {
      chr <- info[2]
      start <- as.numeric(info[3])
      end <- as.numeric(info[4])
    }
    
    return(c(type, chr, start, end))
    
  })
  sv.info.df <- data.frame(t(sv.info.df), stringsAsFactors = FALSE)
  colnames(sv.info.df) <- c("type", "chr", "start", "end")
  
  # create list to store duplicate SVs
  cat("Finding overlapping SVs\n", sep = "")
  duplicates <- list()
  
  for (sv.type in unique(sv.info.df[, "type"])) {
    
    # skip translocations
    if (sv.type != "tra") {
      
      cat("  checking ", sv.type, "\n", sep = "")
      
      # subset by sv type
      sv.info.df.type <- sv.info.df[which(sv.info.df[, "type"] == sv.type), ]
      
      # find SVs from the same type that overlap
      registerDoParallel(num.cores)
      duplicates.by.chr <- foreach(chr=unique(sv.info.df.type[, "chr"]), .combine = c) %dopar% {
        
        # look chromosome by chromosome
        sv.info.df.type.chr <- sv.info.df.type[which(sv.info.df.type[, "chr"] == chr), ]
        # remove svs with exact same id
        sv.info.df.type.chr <- sv.info.df.type.chr[!duplicated(sv.info.df.type.chr), ]
        # sort data by position in chromosome
        sv.info.df.type.chr <- sv.info.df.type.chr[order(as.numeric(sv.info.df.type.chr[, "start"]),
                                                         as.numeric(sv.info.df.type.chr[, "end"])), ]
        # create list to store duplicates of a chromosome
        list.dup.chr <- list()
        
        row = 1
        
        while(row < NROW(sv.info.df.type.chr)) {
          
          # get range of sv
          sv.start <- as.numeric(sv.info.df.type.chr[row, "start"])
          sv.end <- as.numeric(sv.info.df.type.chr[row, "end"])
          # check which svs overlap
          overlap.sv <- which(as.numeric(sv.info.df.type.chr[, "start"]) <= sv.end)
          # make sure to exclude previous svs
          overlap.sv <- overlap.sv[!overlap.sv <= row]
          
          if (length(overlap.sv) > 0) {
            # get duplicate names
            sv.dups <- apply(sv.info.df.type.chr[seq(row, max(overlap.sv), by = 1), ], MARGIN = 1, function(sv) {
              return(paste0(sv, collapse = "."))
            })
            # now add to list
            # duplicates[[length(duplicates)+1]] <- sv.dups
            list.dup.chr[[length(list.dup.chr)+1]] <- sv.dups
            # go to row corresponding to the last duplicate
            row <- max(overlap.sv) + 1
          } else {
            # or go to next line if no duplicates
            row <- row + 1
          }
          
        }
        
        # return list of dups for a chromosome
        list.dup.chr
        
      }
      stopImplicitCluster()
      
      # combine lists from multiple chromosomes
      duplicates <- append(duplicates, duplicates.by.chr)
      
    }
  }
  
  cat("Collapsing and removing SVs\n", sep = "")
  
  # create df to keep collapsed SVs
  sv.parents.collapsed <- data.frame(matrix(nrow = 0, ncol = NCOL(sv.parents)))
  colnames(sv.parents.collapsed) <- colnames(sv.parents)
  
  # create vector to keep SVs to be removed
  svs.to.remove <- c()
  
  # find and collapse duplicates per chromosome
  for (chr in unique(sv.parents[, "chrom"])) {
    
    # subset by chromosome
    sv.parents.chr <- sv.parents[which(sv.parents[, "chrom"] == chr), ]
    
    if (grepl("SCAF", chr) == FALSE) {
      duplicates.chr <- duplicates[grepl(paste0("chr", chr, "."), duplicates, fixed = TRUE)]
    } else {
      duplicates.chr <- duplicates[grepl(paste0(tolower(chr), "."), duplicates, fixed = TRUE)]
    }
    
    if (length(duplicates.chr) > 0) {
      
      cat("  chromosome ", chr, "\n", sep = "")

      # get collapsed rows
      collapsed.rows <- mclapply(duplicates.chr, FUN = function(dup, chr.df) {
        
        dup <- unlist(dup)
        # get all rows that have that position
        chr.df.dup <- chr.df[which(chr.df[, 1] %in% dup), ]
        chr.df.dup <- chr.df.dup[match(dup, chr.df.dup[, 1]), ]
        # create list of SV ids to remove for that chromosome
        svs.to.remove.chr <- c()
        
        for (row in 2:NROW(chr.df.dup)) {
          
          sv1 <- chr.df.dup[row - 1, 12:NCOL(chr.df.dup)]
          sv2 <- chr.df.dup[row, 12:NCOL(chr.df.dup)]
          
          collapsed.svs <- sapply(rbind(sv1, sv2), function(geno) {
            # if there's only one unique genotype called in all duplicates
            if (length(unique(geno)) == 1) {
              # return that genotype
              return(unique(geno))
              # if there's more than one genotype call
            } else {
              # check how many other non-missing genotypes there are
              not.NN.geno <- geno[which(geno != "NN")]
              # if there's only one unique genotype other than NN
              if (length(unique(not.NN.geno)) == 1) {
                # return that genotype
                return(unique(not.NN.geno))
              } else {
                # if there's more than one, it means there's a disagreement between calls,
                # so return "NN"
                return("disagree")
              }
            }
          })
          
          # if there's no disagreement
          if (!"disagree" %in% collapsed.svs) {
            
            # get info about svs
            sv1.info <- chr.df.dup[row - 1, 1:11]
            sv1.id <- unlist(strsplit(chr.df.dup[row - 1, 1], split = ".", fixed = TRUE))
            sv1.size <- as.numeric(sv1.id[4]) - as.numeric(sv1.id[3])
            sv1.missing <- sum(sv1 == "NN")
            
            sv2.info <- chr.df.dup[row, 1:11]
            sv2.id <- unlist(strsplit(chr.df.dup[row, 1], split = ".", fixed = TRUE))
            sv2.size <- as.numeric(sv2.id[4]) - as.numeric(sv2.id[3])
            sv2.missing <- sum(sv2 == "NN")
            
            
            if (sv1.missing < sv2.missing) {
              # if SV1 has less missing data
              sv.to.keep.info <- sv1.info
            } else if (sv1.missing > sv2.missing) {
              # if SV2 has less missing data
              sv.to.keep.info <- sv2.info
            } else {
              # if SVs have the same amount of missing data
              if (sv1.size > sv2.size) {
                # if SV1 is bigger
                sv.to.keep.info <- sv1.info
              } else {
                # if SV2 is bigger
                sv.to.keep.info <- sv2.info
              }
            }
            
            # get the sv id to be removed
            sv.to.remove.id <- c(sv1.info[,1], sv2.info[,1])[c(sv1.info[,1], sv2.info[,1]) != sv.to.keep.info[,1]]
            svs.to.remove.chr <- append(svs.to.remove.chr, sv.to.remove.id)
            
            # make row - 1 NAs
            chr.df.dup[row - 1, ] <- NA
            # make row as collapsed, and get the info for the largest SV
            chr.df.dup[row, ] <- append(sv.to.keep.info, collapsed.svs)
            
          }
          
          # so, if there's a disagreement, keep the way it is
          
        }
        
        # remove rows with NAs
        chr.df.dup <- chr.df.dup[rowSums(is.na(chr.df.dup)) != NCOL(chr.df.dup), ]
        # reorder by middle position
        chr.df.dup <- chr.df.dup[order(chr.df.dup[, "pos"]), ]
        
        # return the df with duplicates after filtering
        return(list(df = chr.df.dup, vector = svs.to.remove.chr))
        
      }, sv.parents.chr, mc.cores = num.cores)
      
      # separate two results from previous loop
      collapsed.rows.df <- do.call(rbind, mclapply(collapsed.rows, function(x) return(x$df), mc.cores = num.cores))
      collapsed.rows.vector <- unlist(lapply(collapsed.rows, function(x) return(x$vector)))
      
      # add results to collapsed df and to the vector with all svs to be removed
      sv.parents.collapsed <- rbind(sv.parents.collapsed, collapsed.rows.df)
      svs.to.remove <- append(svs.to.remove, collapsed.rows.vector)
      
    }
  }

  # now I have to filter the main data.frame (sv.parents) to exclude SVs in the svs.to.remove vector...
  sv.parents <- sv.parents[!sv.parents[, 1] %in% svs.to.remove, ]
  # ...and change the genotype calls for collapsed/kept SVs
  sv.parents <- sv.parents[which(!sv.parents[, 1] %in% sv.parents.collapsed[, 1]), ]
  sv.parents <- rbind(sv.parents, sv.parents.collapsed)
  sv.parents <- SortMixedColumn(sv.parents)

  cat("  ", NROW(sv.parents.collapsed[,1]), " SVs collapsed\n", sep = "")
  cat("  ", length(svs.to.remove), " SVs removed\n", sep = "")
  
  # add one more iteration
  i = i + 1
  
}


# get list of SVs kept
final.svs.to.keep <- sv.parents[, 1]

# load and filter RIL data
cat("Done with iterations, loading RIL data\n", sep = "")
sv.rils <- fread(sv.rils.file, header = TRUE, data.table = FALSE)

cat("Removing duplicated SVs in RIL data\n", sep = "")
sv.rils <- sv.rils[which(sv.rils[, 1] %in% final.svs.to.keep), ]

# write parental and RIL SV filtered files
cat("Writing filtered files\n", sep = "")
outfile.parents <- gsub("hmp.txt", "duplicated-SVs-removed.hmp.txt", sv.parents.file)
outfile.rils <- gsub("hmp.txt", "duplicated-SVs-removed.hmp.txt", sv.rils.file)

fwrite(x = sv.parents, file = outfile.parents, na = NA, quote = FALSE, sep = "\t")
cat("  parents done!\n", sep = "")
fwrite(x = sv.rils, file = outfile.rils, na = NA, quote = FALSE, sep = "\t")
cat("  RILs done!\n", sep = "")


