# # testing
dup.matrix.file <- "tamden_add_private.txt"

# load libraries
library(data.table)

# load matrix
dup.matrix <- fread(dup.matrix.file, header = TRUE, data.table = FALSE, colClasses = 'character')

# identify duplicated IDs for each parent, and update the matrix
for (curr.parent in 1:3) {
  
  # first get a list of IDs
  ID.list <- sapply(dup.matrix[, curr.parent], function(id) {
    id <- unlist(strsplit(id, split = ";"))
    return(id)
  })
  ID.list <- as.character(unlist(ID.list))
  # count number of times each ID appears
  ID.dups <- data.frame(table(ID.list))
  # remove any NAs
  ID.dups <- ID.dups[which(ID.dups[, 1] != "NA"), ]
  # keep only IDs that appear more than once
  ID.dups <- subset(ID.dups, Freq > 1)
  ID.dups <- as.character(ID.dups[,  1])
  
  cat(length(ID.dups), " duplicated IDs in ", colnames(dup.matrix)[curr.parent], "\n", sep = "")
  
  # only proceed if there's duplicated ID for that parent
  if (length(ID.dups) > 0) {
    
    cat("  collapsing...", sep = "")
    
    # for each gene ID that appears more than once
    for (ID in ID.dups) {
      
      rows.to.merge <- grep(ID, dup.matrix[, curr.parent])
      
      merged.ids <- apply(dup.matrix[rows.to.merge, ], MARGIN = 2, function(col) {
        # merge ids in that column
        ids <- c(paste0(col, collapse = ";"))
        # get only unique ids
        ids <- unique(unlist(strsplit(ids, split = ";")))
        # remove NAs
        ids <- ids[ids != "NA"]
        # but if there's no ID left for that column, keep NA
        if (length(ids) == 0) ids <- "NA" 
        # merge unique ids
        ids <- paste0(ids, collapse = ";")
        return(ids)
      })
      
      # remove rows prior merging
      dup.matrix <- dup.matrix[-rows.to.merge, ]
      # add back merged rows at the end of matrix
      dup.matrix <- rbind(dup.matrix, merged.ids)
      
    }
    
    cat(" done!\n", sep = "")
    
  }
}

# write file
outfile.name <- gsub(".csv", ".collapsed.csv", dup.matrix.file)
fwrite(dup.matrix, file = outfile.name, sep = ",", row.names = FALSE, col.names = FALSE, na = "NA", quote = FALSE)
