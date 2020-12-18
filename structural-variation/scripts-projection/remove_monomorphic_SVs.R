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




#### remove monomorphic SVs ----

# load file
sv.parents <- fread(sv.parents.file, header = TRUE, data.table = FALSE)

# get monomorphic SVs that have alternate allele in B73
mono.SVs <- apply(sv.parents, MARGIN = 1, function(sv) {
  
  sv.name <- sv[1]
  geno <- unique(sv[12:length(sv)])
  geno <- geno[geno != "NN"]
  
  if (length(geno) == 1) return(sv.name)
  
})
mono.SVs <- as.character(do.call(c, mono.SVs))

# remove SVs from hmp
sv.parents.no.mono <- sv.parents[which(!sv.parents[, 1] %in% mono.SVs), ]

cat(length(mono.SVs), " monormorphic SVs removed\n", sep = "")

outfile <- gsub("hmp.txt", "poly.hmp.txt", sv.parents.file)
fwrite(x = sv.parents.no.mono, file = outfile, na = NA, quote = FALSE, sep = "\t")

