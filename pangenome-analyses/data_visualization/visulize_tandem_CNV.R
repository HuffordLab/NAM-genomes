library(plotrix)
library(tidyr)
library(matrixStats)
setwd("~/Desktop/pan_genome_nov 2/QC_set/tandem_QC/")

tandem_list <- read.csv("Tandem_list_all_recode_copy_number.txt",sep = "\t",header=FALSE)
colnames(tandem_list) <- c("Query_gene","NAM_genome","Tandem")
head(tandem_list)
#reshape the dataset 

NAM_base <- subset(tandem_list,tandem_list$NAM_genome != "pan_gene" & tandem_list$NAM_genome !="Query_gene" )
#fill_info[rowSums(tandem_list=="")!=ncol(tandem_list), ]
reshaped_datafram <- pivot_wider(NAM_base, names_from = NAM_genome, values_from = Tandem)

# reorder column name before joint the dataset 
matrix_for_count= subset(reshaped_datafram, select = c(Query_gene,B73,Tzi8,Ky21,M162W,Ms71,Oh7B,Oh43,M37W,Mo18W,NC350,HP301,Il14H,P39,CML52,CML69,Ki11,CML228,CML247,CML277,CML322,CML333,Ki3,CML103,Tx303,NC358,B97))
matrix_for_count$tandem_count <- rowSums(matrix_for_count == "Tandem")
table(matrix_for_count$tandem_count)


#tandem copy number matrix 
tandem_copy_list <- read.csv("Tandem_list_all_recode_copy_number.txt",sep = "\t",header=FALSE)
colnames(tandem_copy_list) <- c("Query_gene","NAM_genome","Tandem")
head(tandem_copy_list)
#reshape the dataset 
library(tidyr)
NAM_base_copy_number <- subset(tandem_copy_list,tandem_copy_list$NAM_genome != "pan_gene" & tandem_copy_list$NAM_genome !="Query_gene" )
#fill_info[rowSums(tandem_copy_list=="")!=ncol(tandem_copy_list), ]
reshaped_datafram_copy_number <- pivot_wider(NAM_base_copy_number, names_from = NAM_genome, values_from = Tandem)

# reorder column name before joint the dataset 
matrix_for_copy_count= subset(reshaped_datafram_copy_number, select = c(Query_gene,B73,Tzi8,Ky21,M162W,Ms71,Oh7B,Oh43,M37W,Mo18W,NC350,HP301,Il14H,P39,CML52,CML69,Ki11,CML228,CML247,CML277,CML322,CML333,Ki3,CML103,Tx303,NC358,B97))
# add in pan gene type information, 0 will be considered as NA, missing gene
matrix_for_copy_count$occurance_genome <- rowSums(matrix_for_copy_count[2:27]!=0)



range_matrix <- as.matrix(matrix_for_copy_count[2:27])
copy_number_range <- rowRanges(range_matrix, na.rm = FALSE, dims = 1, n = NULL)
colnames(copy_number_range) <- c("Min_copy","Max_copy")

tanden_gene_type_CNV_range_matrix <- as.data.frame(cbind(matrix_for_copy_count,copy_number_range))

tanden_gene_type_CNV_range_matrix$spread <- tanden_gene_type_CNV_range_matrix$Max_copy-tanden_gene_type_CNV_range_matrix$Min_copy
table(tanden_gene_type_CNV_range_matrix$spread)

write.csv(tanden_gene_type_CNV_range_matrix,file = "tanden_gene_type_CNV_range_matrix.csv")
as.data.frame(table(tanden_gene_type_CNV_range_matrix$spread))
write.csv()
mean(tanden_gene_type_CNV_range_matrix$spread)
# 1.875761
std.error(tanden_gene_type_CNV_range_matrix$spread)

# whole set 
table(tanden_gene_type_CNV_range_matrix$spread)
whole_set <- ggplot(tanden_gene_type_CNV_range_matrix,aes(spread)) + stat_bin(binwidth=1) + #stat_bin(binwidth=1, geom="text", aes(label=..count..), vjust=-1.5) + 
  scale_x_continuous(breaks = seq(0, 135, by = 10)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),text = element_text(size = 20)) + 
  #scale_y_continuous(breaks = seq(0, 10, by = 1))  +
  ylab("Number of Putative Tandem Duplicates") + xlab("Putative Tandem Duplicates Copy Number Variation (0-135)") + theme_classic()


#subset for zoom in 
copy_variation_10_and_below <- subset(tanden_gene_type_CNV_range_matrix,tanden_gene_type_CNV_range_matrix$spread <=10)
table(copy_variation_10_and_below$spread )
# creating grobs to add in 
bp<- ggplotGrob(ggplot(copy_variation_10_and_below,aes(spread)) + stat_bin(binwidth=1) + #stat_bin(binwidth=1, geom="text", aes(label=..count..), vjust=-0.5,angle =0,colour="red",size=3) + 
                  scale_x_continuous(breaks = seq(0, 10, by = 1)) + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),text = element_text(size = 12)) + theme_classic() +
                  ylab("Number of Putative Tandem Duplicates") +  xlab("Putative Tandem Duplicates Copy Number Variation (0-10)"))


tanden_cnv <- whole_set + 
  annotation_custom(
    grob = bp,
    xmin = 35,
    xmax = 125,
    ymin = 1000,
    ymax = 7500
  ) 
tanden_cnv + theme(text = element_text(size = 14))


# count tandem specific to a genome
tandem_type <- read.csv("tandem_type_stats.csv")
table(tandem_type$Tandem)
#   1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20 
#6556 2743 1411  925  556  424  348  274  205  214  147  119  125  108   85   73   79   80   65   77 
#  21   22   23   24   25   26 
#  80   88  119  125  271  970 
