# bootstraping 
march_pan_matrix <- read.csv(file = "pan_gene_matrix_march_v3_all_info.csv",header = TRUE,stringsAsFactors=FALSE)

# replace the gene ID using pan gene ID to speed up the bootstrapping 
pan_gene_id_df <- sapply(march_pan_matrix[,-1], function(x) {ind <- which(x!="NA"); x[ind] = march_pan_matrix[ind,1]; return(x)})


# bootstrap for all 26 genomes 
all_26_genomes <- pan_gene_id_df[,3:28]
col_names <- colnames(all_26_genomes)
for (i in 1 : ncol(all_26_genomes)) {
  nam_pan_gene_id <- as.matrix(all_26_genomes[,i])
  write.csv(nam_pan_gene_id, file = paste(col_names[i], "_pan_id.csv", sep = ""))
}

# bootstrap for B73 B97 Ky21 M162W Ms71 Oh43 and Oh7B
temperate = as.data.frame(all_26_genomes) %>% select("B73","B97", "Ky21", "M162W", "Ms71", "Oh43","Oh7B","HP301","P39","Il14H")
# remove pan genes that are present in non of the 7 genomes 
temperate$genome_presence <- 10 - apply(temperate, 1, function(x) sum(is.na(x)))
temperate_pan_matrix =temperate %>% filter(genome_presence >0) 
# there is 81,097 pan genes for these 10 lines 

col_names <- colnames(temperate_pan_matrix)
for (i in 1 : ncol(temperate_pan_matrix)) {
  nam_pan_gene_id <- as.matrix(temperate_pan_matrix[,i])
  write.csv(nam_pan_gene_id, file = paste(col_names[i], "_temperate_pan_id.csv", sep = ""))
}

# bootstrap for B73 CML52 CML69 CML103 CML228 CML247 CML277 CML322 CML333 Ki3 Ki11 NC350 NC358, and Tzi8
tropical = as.data.frame(all_26_genomes) %>% select("CML52", "CML69", "CML103", "CML228", "CML247","CML277","CML322","CML333","Ki3","Ki11","NC350","NC358","Tzi8")
# remove pan genes that are present in non of the 7 genomes 
tropical$genome_presence <- 13 - apply(tropical, 1, function(x) sum(is.na(x)))
tropical_pan_matrix =tropical %>% filter(genome_presence >0) 
# there is 85923 pan genes for these 13 lines 

col_names <- colnames(tropical_pan_matrix)
for (i in 1 : ncol(tropical_pan_matrix)) {
  nam_pan_gene_id <- as.matrix(tropical_pan_matrix[,i])
  write.csv(nam_pan_gene_id, file = paste(col_names[i], "_tropical_pan_id.csv", sep = ""))
}

