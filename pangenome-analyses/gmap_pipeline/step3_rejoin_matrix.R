setwd("~/Desktop/pan_genome_nov/")
fill_info <- read.csv("CDS_canonical.txt",sep = "\t",header=FALSE)
colnames(fill_info) <- c("Query_gene","NAM_genome","Gmap_coordinate")

#reshape the dataset 
library(tidyr)
fill_info[rowSums(fill_info=="")!=ncol(fill_info), ]
reshaped_datafram <- pivot_wider(fill_info, names_from = NAM_genome, values_from = Gmap_coordinate)
write.csv(reshaped_datafram, file= "reshaped_matrix_final.csv")

# reorder column name before joint the dataset 
matrix_for_merge = subset(reshaped_datafram, select = c(Query_gene,B73,Tzi8,Ky21,M162W,Ms71,Oh7B,Oh43,M37W,Mo18W,NC350,HP301,Il14H,P39,CML52,CML69,Ki11,CML228,CML247,CML277,CML322,CML333,Ki3,CML103,Tx303,NC358,B97))
write.csv(matrix_for_merge,file = "~/Desktop/pan_genome_nov 2/ordered_NA_fill_matrix.csv")

# read existing matrix 
pan_26_matrix <- read.csv("pan26_all.collapsed.csv")
pan_26_for_merge = subset(pan_26_matrix, select = c(Query_gene,B73,Tzi8,Ky21,M162W,Ms71,Oh7B,Oh43,M37W,Mo18W,NC350,HP301,Il14H,P39,CML52,CML69,Ki11,CML228,CML247,CML277,CML322,CML333,Ki3,CML103,Tx303,NC358,B97))

merged <- rbind(pan_26_for_merge,matrix_for_merge)
write.csv(merged,file="gmap_join_matrix_final.csv")
# load libraries
library(data.table)


dup.matrix.file <- "gmap_join_matrix_final.csv"

# load libraries
library(data.table)

# load matrix
dup.matrix <- fread(dup.matrix.file, header = TRUE, data.table = FALSE, colClasses = 'character')
# load matrix

# identify duplicated IDs for each parent, and update the matrix
for (curr.parent in 1:2) {
  
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
outfile.name <- gsub(".csv", ".gmap_no_compression_final.csv", dup.matrix.file)
fwrite(dup.matrix, file = outfile.name, sep = ",", row.names = FALSE, col.names = FALSE, na = "NA", quote = FALSE)
