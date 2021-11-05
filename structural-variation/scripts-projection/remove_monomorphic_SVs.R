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
# you should provide 1 argument
if (length(args) != 1) {
  stop("incorrect number of arguments provided.

       Usage: Rscript remove_monomorphic_SVs.R [parents_file]
       ")
}

# assign arguments to variables
sv.parents.file <- args[1]

# sv.parents.file <- "data/NAM_founders_SVs.not-collapsed.hmp.txt"


#### libraries ----

if(!require("data.table")) install.packages("data.table")


removeMonoSVs <- function(hmp, sv_to_remove_info) {
  
  hmp.filtered <- data.frame(stringsAsFactors = FALSE)
  for (chr in unique(hmp[, "chrom"])) {
    
    # subset by chromosome
    hmp.chr <- subset(hmp, chrom == chr)
    sv.info.chr <- subset(sv_to_remove_info, chrom == chr)
    
    rows.to.exclude <- c()
    for (row in 1:NROW(sv.info.chr)) {
      
      # some SV names (especially translocations) may be duplicated,
      # so I had to filter by both name and positions
      rows.to.exclude <- append(rows.to.exclude, which(hmp.chr[, 1] == sv.info.chr[row, 1] & hmp.chr[, 4] == sv.info.chr[row, 3]))
      
    }
    
    # filter hapmap and add to new df
    rows.to.exclude <- sort(unique(rows.to.exclude))
    hmp.filtered <- rbind(hmp.filtered, hmp.chr[-rows.to.exclude, ])
    
  }
  
  return(hmp.filtered)
  
}




#### remove monomorphic SVs ----

# load file
sv.parents <- fread(sv.parents.file, header = TRUE, data.table = FALSE)

# remove SVs in which B73 is TT
sv.parents.no.TT <- sv.parents[which(sv.parents[, "B73"] != "TT"), ]

cat(NROW(sv.parents) - NROW(sv.parents.no.TT), " SVs removed for being called present in B73\n", sep = "")

# get monomorphic SVs
mono.SVs <- apply(sv.parents.no.TT, MARGIN = 1, function(sv) {

  sv.info <- sv[c(1, 3, 4)]
  geno <- unique(sv[12:length(sv)])
  geno <- geno[geno != "NN"]
  
  if (length(geno) <= 1) return(sv.info)

})

# remove extra whitespace on any cell
mono.SVs <- lapply(mono.SVs, function(x) gsub(" +", "", x, perl = TRUE))

# make it data frame
mono.SVs <- data.frame(do.call(rbind, mono.SVs), stringsAsFactors = FALSE)

# remove mono SVs from hmp
sv.parents.no.mono <- removeMonoSVs(sv.parents.no.TT, mono.SVs)

cat(NROW(sv.parents.no.TT) - NROW(sv.parents.no.mono), " monormorphic SVs removed\n", sep = "")

outfile <- gsub("hmp.txt", "poly.hmp.txt", sv.parents.file)
fwrite(x = sv.parents.no.mono, file = outfile, na = NA, quote = FALSE, sep = "\t")
