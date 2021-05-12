library(data.table)

# set up anchor for pan-gene type 
pan_gene_type <- read.csv(file = "pan_gene_matrix_march_v3_all_info.csv")
pan_gene_type_anchor <- pan_gene_type[,c(1,31)]
names(pan_gene_type_anchor)[1] <- "PanGeneID"


files <- list.files(path="NAM_pangene_expression_counts_per_tissue-RPKM/", pattern="*.tsv", full.names=TRUE, recursive=FALSE)
lapply(files, function(x) {
  df <- fread(x, header = TRUE, data.table = FALSE)
  # get sample names based on column names
  samples <- sapply(colnames(df)[8:NCOL(df)], function(name) {
    rep_name <- unlist(strsplit(name, split = "_"))
    rep_name <- rep_name[-length(rep_name)]
    rep_name <- paste0(rep_name, collapse = "_")
    return(rep_name)
  })
  samples <- unique(samples)
  # calculate average for each tissue type
  datalist = list()
  for (name in samples) {
    reps <- data.frame(df[, grep(name, colnames(df))])
    #reps lapply(stam[,4:15], function(x) as.numeric(as.character(x)))
    sample_mean <- round(rowMeans(reps), digits = 4)
    #df_sample <- cbind(df[, 1:2], sample_mean)
    df_name <- as.data.frame(sample_mean)
    colnames(df_name)[1] <- name
    datalist[[name]] <- df_name
  }
  # reshape the list into matrix 
  df_ave <- data.frame(cbind(df[, 1:2], as.data.frame(datalist)))
  df_ave_sub = subset(df_ave, select = - c(grep("embryo", names(df_ave)),grep("endosperm", names(df_ave))))
  df_ave_sub[df_ave_sub < 1] <- NA
  df_ave_sub$number_of_tissue_expression <- 8 - apply(df_ave_sub, 1, function(x) sum(is.na(x)))
  left_join(df_ave_sub,pan_gene_type_anchor) %>% write.csv( paste0(x, ".csv"), row.names = FALSE)
})


# because CML52 and CML277 has less tissue type, they are processed separately with modified script 

# CML277 
# get mean value for expression for each tissue 
# load table
df <- fread("/NAM_pangene_expression_counts_per_tissue-RPKM/CML52_CML277/CML277_full-rpkm.tsv", header = TRUE, data.table = FALSE)
# get sample names based on column names
samples <- sapply(colnames(df)[8:NCOL(df)], function(name) {
  rep_name <- unlist(strsplit(name, split = "_"))
  rep_name <- rep_name[-length(rep_name)]
  rep_name <- paste0(rep_name, collapse = "_")
  return(rep_name)
})
samples <- unique(samples)
# calculate average for each tissue type
datalist = list()
for (name in samples) {
  reps <- data.frame(df[, grep(name, colnames(df))])
  sample_mean <- round(rowMeans(reps), digits = 4)
  #df_sample <- cbind(df[, 1:2], sample_mean)
  df_name <- as.data.frame(sample_mean)
  colnames(df_name)[1] <- name
  datalist[[name]] <- df_name
}
# reshape the list into matrix 
RPKM_CML277 <- cbind(df[, 1:2], as.data.frame(datalist))
RPKM_CML277_sub = RPKM_CML277
RPKM_CML277_sub[RPKM_CML277_sub < 1] <- NA
RPKM_CML277_sub$number_of_tissue_expression <- 8 - apply(RPKM_CML277_sub, 1, function(x) sum(is.na(x)))
CML277_Tissue_type_count <- left_join(RPKM_CML277_sub,pan_gene_type_anchor) %>% write.csv("CML277_full-rpkm.tsv.csv")



# looping through all the files in the folder 
# CML52 
# get mean value for expression for each tissue 
# load table
df <- fread("NAM_pangene_expression_counts_per_tissue-RPKM/CML52_CML277/CML52_full-RPKM.tsv", header = TRUE, data.table = FALSE)
# get sample names based on column names
samples <- sapply(colnames(df)[8:NCOL(df)], function(name) {
  rep_name <- unlist(strsplit(name, split = "_"))
  rep_name <- rep_name[-length(rep_name)]
  rep_name <- paste0(rep_name, collapse = "_")
  return(rep_name)
})
samples <- unique(samples)
# calculate average for each tissue type
datalist = list()
for (name in samples) {
  reps <- data.frame(df[, grep(name, colnames(df))])
  sample_mean <- round(rowMeans(reps), digits = 4)
  #df_sample <- cbind(df[, 1:2], sample_mean)
  df_name <- as.data.frame(sample_mean)
  colnames(df_name)[1] <- name
  datalist[[name]] <- df_name
}
# reshape the list into matrix 
RPKM_CML52 <- cbind(df[, 1:2], as.data.frame(datalist))
RPKM_CML52_sub = RPKM_CML52
RPKM_CML52_sub[RPKM_CML52_sub < 1] <- NA
RPKM_CML52_sub$number_of_tissue_expression <- 8 - apply(RPKM_CML52_sub, 1, function(x) sum(is.na(x)))
CML52_Tissue_type_count <- left_join(RPKM_CML52_sub,pan_gene_type_anchor) 
write.csv(CML52_Tissue_type_count,"CML52_full-RPKM.tsv.csv")

