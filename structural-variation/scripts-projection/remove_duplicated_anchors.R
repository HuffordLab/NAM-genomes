library(data.table)

# assign arguments to variables
args <- commandArgs(trailingOnly = TRUE)
hmp_file <- args[1]

# remove duplicates
hmp <- fread(hmp_file, header = TRUE, data.table = FALSE)
hmp <- hmp[!duplicated(hmp[, 1]) & !duplicated(hmp[, 4]), ]

# write filtered file
fwrite(hmp, hmp_file, quote = FALSE, sep = "\t", na = NA, row.names = FALSE)
